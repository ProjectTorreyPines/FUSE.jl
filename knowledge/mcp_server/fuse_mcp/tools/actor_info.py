"""MCP tools for querying FUSE actor information"""

from typing import Any, Dict, List, Optional
from mcp.types import Tool
from mcp.types import TextContent
import json

from ..actor_knowledge import ActorKnowledgeBase


class ActorInfoTools:
    """Collection of MCP tools for actor information queries"""
    
    def __init__(self, knowledge_base: ActorKnowledgeBase):
        self.kb = knowledge_base
    
    def get_tools(self) -> List[Tool]:
        """Return all actor information tools"""
        return [
            Tool(
                name="list_actors",
                description="List all available FUSE actors, optionally filtered by category",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "category": {
                            "type": "string",
                            "description": "Optional category to filter by (e.g., 'transport', 'equilibrium', 'hcd')"
                        }
                    }
                }
            ),
            Tool(
                name="get_actor_info",
                description="Get detailed information about a specific FUSE actor",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "actor_name": {
                            "type": "string",
                            "description": "Name of the actor to get information about"
                        }
                    },
                    "required": ["actor_name"]
                }
            ),
            Tool(
                name="search_actors",
                description="Search for actors by name, description, or category",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "query": {
                            "type": "string",
                            "description": "Search query to find relevant actors"
                        }
                    },
                    "required": ["query"]
                }
            ),
            Tool(
                name="get_actor_categories",
                description="List all available actor categories",
                inputSchema={
                    "type": "object",
                    "properties": {}
                }
            ),
            Tool(
                name="get_actor_usage_guide",
                description="Get a comprehensive usage guide for a specific actor",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "actor_name": {
                            "type": "string",
                            "description": "Name of the actor to get usage guide for"
                        }
                    },
                    "required": ["actor_name"]
                }
            )
        ]
    
    async def handle_tool_call(self, name: str, arguments: Dict[str, Any]) -> List[TextContent]:
        """Handle tool calls for actor information"""
        
        if name == "list_actors":
            category = arguments.get("category")
            actors = self.kb.list_actors(category)
            
            if category:
                response = f"Actors in category '{category}':\n"
            else:
                response = "All available actors:\n"
            
            response += "\n".join(f"- {actor}" for actor in sorted(actors))
            return [TextContent(type="text", text=response)]
        
        elif name == "get_actor_info":
            actor_name = arguments["actor_name"]
            actor = self.kb.get_actor(actor_name)
            
            if not actor:
                return [TextContent(type="text", text=f"Actor '{actor_name}' not found.")]
            
            response = f"# {actor.name} Actor\n\n"
            response += f"**Category**: {actor.category}\n"
            response += f"**Description**: {actor.description}\n\n"
            
            if actor.parameters:
                response += "## Parameters\n"
                for param in actor.parameters:
                    response += f"- **{param.name}** ({param.type}, {param.unit}): {param.description}\n"
                    if param.default is not None:
                        response += f"  Default: {param.default}\n"
                response += "\n"
            
            if actor.inputs:
                response += f"**Inputs**: {', '.join(actor.inputs)}\n"
            
            if actor.outputs:
                response += f"**Outputs**: {', '.join(actor.outputs)}\n"
            
            return [TextContent(type="text", text=response)]
        
        elif name == "search_actors":
            query = arguments["query"]
            results = self.kb.search_actors(query)
            
            if not results:
                return [TextContent(type="text", text=f"No actors found matching '{query}'.")]
            
            response = f"Actors matching '{query}':\n\n"
            for actor in results:
                response += f"**{actor.name}** ({actor.category}): {actor.description}\n"
            
            return [TextContent(type="text", text=response)]
        
        elif name == "get_actor_categories":
            categories = self.kb.get_categories()
            response = "Available actor categories:\n"
            response += "\n".join(f"- {cat}" for cat in sorted(categories))
            return [TextContent(type="text", text=response)]
        
        elif name == "get_actor_usage_guide":
            actor_name = arguments["actor_name"]
            guide = self.kb.get_actor_usage_guide(actor_name)
            return [TextContent(type="text", text=guide)]
        
        else:
            return [TextContent(type="text", text=f"Unknown tool: {name}")]