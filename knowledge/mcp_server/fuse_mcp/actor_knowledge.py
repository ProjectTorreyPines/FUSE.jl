"""Actor knowledge management system for FUSE MCP server"""

import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from loguru import logger


@dataclass
class ActorParameter:
    """Represents a single actor parameter"""
    name: str
    type: str
    unit: str
    description: str
    default: Optional[Any] = None
    required: bool = True


@dataclass
class ActorInfo:
    """Complete information about a FUSE actor"""
    name: str
    category: str
    description: str
    parameters: List[ActorParameter]
    inputs: List[str]
    outputs: List[str]
    file_path: str
    julia_type: Optional[str] = None
    usage_examples: Optional[List[str]] = None


class ActorKnowledgeBase:
    """Manages FUSE actor knowledge from JSON files"""
    
    def __init__(self, actors_dir: Path, knowledge_base_file: Path):
        self.actors_dir = Path(actors_dir)
        self.knowledge_base_file = Path(knowledge_base_file)
        self.actors: Dict[str, ActorInfo] = {}
        self.categories: Dict[str, List[str]] = {}
        self._load_actors()
    
    def _load_actors(self):
        """Load all actor definitions from JSON files"""
        logger.info(f"Loading actors from {self.actors_dir}")

        # Load from individual JSON files if they exist
        if self.actors_dir.exists():
            for json_file in self.actors_dir.rglob("*.json"):
                try:
                    with open(json_file, 'r') as f:
                        actor_data = json.load(f)

                    actor_info = self._parse_actor_json(actor_data, json_file)
                    if actor_info:
                        self.actors[actor_info.name] = actor_info

                        # Organize by category
                        if actor_info.category not in self.categories:
                            self.categories[actor_info.category] = []
                        self.categories[actor_info.category].append(actor_info.name)

                except Exception as e:
                    logger.warning(f"Failed to load actor from {json_file}: {e}")

        # Load from unified knowledge base if it exists
        if self.knowledge_base_file.exists():
            self._load_unified_knowledge_base()

        logger.info(f"Loaded {len(self.actors)} actors across {len(self.categories)} categories")
    
    def _parse_actor_json(self, data: Dict, file_path: Path) -> Optional[ActorInfo]:
        """Parse actor JSON data into ActorInfo object"""
        try:
            # Extract basic info
            name = data.get("name", file_path.stem.replace("_actor", ""))
            description = data.get("description", "No description available")
            
            # Determine category from file path
            category_parts = file_path.parent.name
            if category_parts == "actors":
                category = "general"
            else:
                category = category_parts
            
            # Parse parameters - handle both formats
            parameters = []
            params_data = data.get("parameters", {}) or data.get("key_parameters", {})
            for param_name, param_info in params_data.items():
                if isinstance(param_info, dict):
                    param = ActorParameter(
                        name=param_name,
                        type=param_info.get("type", "unknown"),
                        unit=param_info.get("unit", "-"),
                        description=param_info.get("description", "No description"),
                        default=param_info.get("default"),
                        required=param_info.get("required", True)
                    )
                    parameters.append(param)
                elif isinstance(param_info, str):
                    # Simple format where value is just description
                    param = ActorParameter(
                        name=param_name,
                        type="unknown",
                        unit="-",
                        description=param_info,
                        default=None,
                        required=True
                    )
                    parameters.append(param)
            
            # Extract inputs/outputs - handle different field names
            inputs = data.get("inputs", []) or data.get("data_inputs", [])
            outputs = data.get("outputs", []) or data.get("data_outputs", [])
            
            return ActorInfo(
                name=name,
                category=category,
                description=description,
                parameters=parameters,
                inputs=inputs,
                outputs=outputs,
                file_path=str(file_path),
                julia_type=data.get("julia_type"),
                usage_examples=data.get("usage_examples")
            )
        
        except Exception as e:
            logger.error(f"Error parsing actor JSON from {file_path}: {e}")
            return None
    
    def _load_unified_knowledge_base(self):
        """Load additional data from unified knowledge base"""
        try:
            with open(self.knowledge_base_file, 'r') as f:
                kb_data = json.load(f)

            actors_data = kb_data.get("actors", {})

            # Load actors from unified knowledge base
            for actor_name, actor_info in actors_data.items():
                if actor_name not in self.actors:
                    # Create new actor from knowledge base
                    parameters = []
                    # Handle both 'parameters' and 'key_parameters' keys
                    params_data = actor_info.get("parameters", actor_info.get("key_parameters", {}))
                    for param_name, param_data in params_data.items():
                        if isinstance(param_data, dict):
                            param = ActorParameter(
                                name=param_name,
                                type=param_data.get("type", "unknown"),
                                unit=param_data.get("unit", "-"),
                                description=param_data.get("description", "No description"),
                                default=param_data.get("default"),
                                required=param_data.get("required", True)
                            )
                            parameters.append(param)
                        elif isinstance(param_data, str):
                            # Simple string description
                            param = ActorParameter(
                                name=param_name,
                                type="unknown",
                                unit="-",
                                description=param_data,
                                default=None,
                                required=True
                            )
                            parameters.append(param)

                    category = actor_info.get("category", "general")
                    # Handle both 'inputs'/'outputs' and 'data_inputs'/'data_outputs'
                    inputs = actor_info.get("inputs", actor_info.get("data_inputs", []))
                    outputs = actor_info.get("outputs", actor_info.get("data_outputs", []))

                    new_actor = ActorInfo(
                        name=actor_name,
                        category=category,
                        description=actor_info.get("description", "No description available"),
                        parameters=parameters,
                        inputs=inputs,
                        outputs=outputs,
                        file_path=str(self.knowledge_base_file),
                        julia_type=actor_info.get("julia_type"),
                        usage_examples=actor_info.get("usage_examples")
                    )

                    self.actors[actor_name] = new_actor

                    # Organize by category
                    if category not in self.categories:
                        self.categories[category] = []
                    self.categories[category].append(actor_name)
                else:
                    # Enhance existing actors with additional information
                    if "usage_examples" in actor_info:
                        self.actors[actor_name].usage_examples = actor_info["usage_examples"]

        except Exception as e:
            logger.warning(f"Failed to load unified knowledge base: {e}")
    
    def get_actor(self, name: str) -> Optional[ActorInfo]:
        """Get actor info by name"""
        return self.actors.get(name)
    
    def list_actors(self, category: Optional[str] = None) -> List[str]:
        """List all actor names, optionally filtered by category"""
        if category:
            return self.categories.get(category, [])
        return list(self.actors.keys())
    
    def get_categories(self) -> List[str]:
        """Get all available actor categories"""
        return list(self.categories.keys())
    
    def search_actors(self, query: str) -> List[ActorInfo]:
        """Search actors by name or description"""
        query_lower = query.lower()
        results = []
        
        for actor in self.actors.values():
            if (query_lower in actor.name.lower() or 
                query_lower in actor.description.lower() or
                query_lower in actor.category.lower()):
                results.append(actor)
        
        return results
    
    def get_actor_usage_guide(self, actor_name: str) -> str:
        """Generate a usage guide for a specific actor"""
        actor = self.get_actor(actor_name)
        if not actor:
            return f"Actor '{actor_name}' not found."
        
        guide = f"# {actor.name} Actor Usage Guide\n\n"
        guide += f"**Category**: {actor.category}\n"
        guide += f"**Description**: {actor.description}\n\n"
        
        if actor.parameters:
            guide += "## Parameters\n\n"
            for param in actor.parameters:
                guide += f"- **{param.name}** ({param.type}): {param.description}\n"
                guide += f"  - Unit: {param.unit}\n"
                if param.default is not None:
                    guide += f"  - Default: {param.default}\n"
                guide += f"  - Required: {'Yes' if param.required else 'No'}\n\n"
        
        if actor.inputs:
            guide += "## Inputs\n"
            for inp in actor.inputs:
                guide += f"- {inp}\n"
            guide += "\n"
        
        if actor.outputs:
            guide += "## Outputs\n"
            for out in actor.outputs:
                guide += f"- {out}\n"
            guide += "\n"
        
        if actor.usage_examples:
            guide += "## Usage Examples\n"
            for example in actor.usage_examples:
                guide += f"```julia\n{example}\n```\n\n"
        
        return guide