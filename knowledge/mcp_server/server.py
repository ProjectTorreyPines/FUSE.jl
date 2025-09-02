#!/usr/bin/env python3
"""
FUSE MCP Server - Model Context Protocol server for FUSE plasma simulation actors

This server provides LLM access to FUSE actor documentation, usage guides, parameter information,
and relationships between actors. It supports comprehensive queries about plasma simulation workflows.
"""

import asyncio
import sys
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

from mcp.server import Server, InitializationOptions
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent, Resource
from fuse_mcp.actor_knowledge import ActorKnowledgeBase


# Global knowledge base
kb = None
relationships_data = None

def load_relationships():
    """Load relationships from the unified knowledge base"""
    global relationships_data
    import json
    
    server_dir = Path(__file__).parent
    knowledge_file = server_dir.parent / "fuse_knowledge_base.json"
    
    try:
        with open(knowledge_file, 'r') as f:
            data = json.load(f)
            relationships_data = data.get("relationships", {})
            print(f"Loaded relationships: {list(relationships_data.keys())}")
    except Exception as e:
        print(f"Failed to load relationships: {e}")
        relationships_data = {}

async def handle_relationships_tool(arguments: dict):
    """Handle the get_actor_relationships tool"""
    if relationships_data is None:
        return [TextContent(type="text", text="❌ Relationships data not loaded")]
    
    actor_name = arguments.get("actor_name")
    relationship_type = arguments.get("relationship_type", "all")
    
    # If no specific actor, show overview
    if not actor_name:
        return [TextContent(type="text", text=generate_relationships_overview())]
    
    # Check if actor exists
    if actor_name not in kb.actors:
        return [TextContent(type="text", text=f"❌ Actor '{actor_name}' not found. Try using 'search_actors' first.")]
    
    # Generate relationships report for specific actor
    response = f"# {actor_name} Relationships\n\n"
    
    if relationship_type in ["data_flow", "all"]:
        response += generate_data_flow_section(actor_name)
    
    if relationship_type in ["dependencies", "all"]:
        response += generate_dependencies_section(actor_name)
    
    if relationship_type in ["physics_groups", "all"]:
        response += generate_physics_groups_section(actor_name)
    
    if relationship_type in ["category_connections", "all"]:
        response += generate_category_connections_section(actor_name)
    
    return [TextContent(type="text", text=response)]

def generate_relationships_overview():
    """Generate overview of all relationship types"""
    response = "# 🔗 FUSE Actor Relationships Overview\n\n"
    
    if "data_flow_connections" in relationships_data:
        data_flows = relationships_data["data_flow_connections"]
        response += f"## 📊 Data Flow Connections\n"
        response += f"- **{len(data_flows)}** actors with data flow information\n"
        
        # Count total connections
        total_connections = 0
        for actor_connections in data_flows.values():
            total_connections += len(actor_connections.get("feeds_to", []))
            total_connections += len(actor_connections.get("receives_from", []))
        
        response += f"- **{total_connections}** total data flow connections\n\n"
    
    if "physics_domain_groups" in relationships_data:
        physics_groups = relationships_data["physics_domain_groups"]
        response += f"## ⚛️ Physics Domain Groups\n"
        for group_name, actors in physics_groups.items():
            response += f"- **{group_name}**: {len(actors)} actors\n"
        response += "\n"
    
    if "sub_actor_dependencies" in relationships_data:
        dependencies = relationships_data["sub_actor_dependencies"]
        response += f"## 🔧 Actor Dependencies\n"
        response += f"- **{len(dependencies)}** actors with sub-actor dependencies\n\n"
    
    if "category_connections" in relationships_data:
        cat_connections = relationships_data["category_connections"]
        response += f"## 📁 Category Connections\n"
        for category, connections in cat_connections.items():
            total_connections = sum(connections.values()) if isinstance(connections, dict) else 0
            response += f"- **{category}**: {total_connections} total connections to {len(connections)} other categories\n"
        response += "\n"
    
    response += "💡 **Usage**: Use `get_actor_relationships` with a specific actor name to see detailed connections."
    
    return response

def generate_data_flow_section(actor_name: str):
    """Generate data flow section for specific actor"""
    response = "## 📊 Data Flow Connections\n\n"
    
    data_flows = relationships_data.get("data_flow_connections", {})
    actor_flows = data_flows.get(actor_name, {})
    
    # What this actor feeds to
    feeds_to = actor_flows.get("feeds_to", [])
    if feeds_to:
        response += "### ➡️ Feeds Data To:\n"
        for connection in feeds_to:
            target_actor = connection.get("actor", "Unknown")
            shared_data = connection.get("shared_data", [])
            response += f"- **{target_actor}**\n"
            for data in shared_data:
                data_path = data.get("data_path", "unknown")
                response += f"  - `{data_path}`\n"
        response += "\n"
    else:
        response += "### ➡️ Feeds Data To: None\n\n"
    
    # What feeds to this actor
    receives_from = actor_flows.get("receives_from", [])
    if receives_from:
        response += "### ⬅️ Receives Data From:\n"
        for connection in receives_from:
            source_actor = connection.get("actor", "Unknown")
            shared_data = connection.get("shared_data", [])
            response += f"- **{source_actor}**\n"
            for data in shared_data:
                data_path = data.get("data_path", "unknown")
                response += f"  - `{data_path}`\n"
        response += "\n"
    else:
        response += "### ⬅️ Receives Data From: None\n\n"
    
    return response

def generate_dependencies_section(actor_name: str):
    """Generate dependencies section for specific actor"""
    response = "## 🔧 Sub-Actor Dependencies\n\n"
    
    dependencies = relationships_data.get("sub_actor_dependencies", {})
    actor_deps = dependencies.get(actor_name, {})
    
    if actor_deps:
        response += f"**{actor_name}** depends on the following sub-actors:\n"
        for dep_type, dep_actors in actor_deps.items():
            response += f"- **{dep_type}**: {', '.join(dep_actors)}\n"
        response += "\n"
    else:
        response += f"**{actor_name}** has no sub-actor dependencies (it's a simple actor).\n\n"
    
    return response

def generate_physics_groups_section(actor_name: str):
    """Generate physics groups section for specific actor"""
    response = "## ⚛️ Physics Domain Groups\n\n"
    
    physics_groups = relationships_data.get("physics_domain_groups", {})
    
    actor_groups = []
    for group_name, actors in physics_groups.items():
        if actor_name in actors:
            actor_groups.append((group_name, actors))
    
    if actor_groups:
        for group_name, group_actors in actor_groups:
            response += f"**{actor_name}** belongs to the **{group_name}** physics domain group:\n"
            other_actors = [a for a in group_actors if a != actor_name]
            if other_actors:
                response += f"- Related actors: {', '.join(other_actors)}\n"
            response += "\n"
    else:
        response += f"**{actor_name}** doesn't belong to any specific physics domain group.\n\n"
    
    return response

def generate_category_connections_section(actor_name: str):
    """Generate category connections section for specific actor"""
    response = "## 📁 Category Connections\n\n"
    
    # Get actor's category
    actor_info = kb.get_actor(actor_name)
    if not actor_info:
        response += "Unable to determine actor category.\n\n"
        return response
    
    actor_category = actor_info.category
    cat_connections = relationships_data.get("category_connections", {})
    
    if actor_category in cat_connections:
        connections = cat_connections[actor_category]
        response += f"**{actor_name}** is in the **{actor_category}** category, which connects to:\n"
        for target_category, connection_count in connections.items():
            response += f"- **{target_category}** ({connection_count} connections)\n"
        response += "\n"
    else:
        response += f"No category connection data available for **{actor_category}**.\n\n"
    
    return response

async def main():
    global kb
    
    # Initialize knowledge base
    server_dir = Path(__file__).parent
    actors_dir = server_dir.parent / "actors"
    knowledge_file = server_dir.parent / "fuse_knowledge_base.json"
    
    print(f"Loading actors from: {actors_dir}")
    print(f"Knowledge base file: {knowledge_file}")
    
    kb = ActorKnowledgeBase(
        actors_dir=actors_dir,
        knowledge_base_file=knowledge_file
    )
    
    # Load relationships data
    load_relationships()
    
    # Create server
    server = Server(name="fuse-mcp", version="0.1.0")
    
    @server.list_tools()
    async def list_tools():
        return [
            Tool(
                name="list_actors",
                description="List all available FUSE actors, optionally filtered by category",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "category": {
                            "type": "string",
                            "description": "Optional category to filter by (e.g., 'transport', 'equilibrium', 'hcd', 'pedestal')"
                        }
                    }
                }
            ),
            Tool(
                name="get_actor_info",
                description="Get detailed information about a specific FUSE actor including parameters, inputs, and outputs",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "actor_name": {
                            "type": "string",
                            "description": "Name of the actor to get information about (e.g., 'ActorTGLF', 'ActorDynamicPlasma')"
                        }
                    },
                    "required": ["actor_name"]
                }
            ),
            Tool(
                name="search_actors",
                description="Search for FUSE actors by name, description, or physics domain",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "query": {
                            "type": "string",
                            "description": "Search query to find relevant actors (e.g., 'transport', 'turbulence', 'TGLF')"
                        }
                    },
                    "required": ["query"]
                }
            ),
            Tool(
                name="get_actor_categories",
                description="List all available FUSE actor categories and physics domains",
                inputSchema={
                    "type": "object",
                    "properties": {}
                }
            ),
            Tool(
                name="get_actor_usage_guide",
                description="Get a comprehensive usage guide for a specific FUSE actor with parameters and examples",
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
            ),
            Tool(
                name="get_actor_relationships",
                description="Get relationships and connections between FUSE actors including data flow, dependencies, and physics domain groupings",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "actor_name": {
                            "type": "string", 
                            "description": "Name of the actor to get relationships for (optional - if not provided, shows overview)"
                        },
                        "relationship_type": {
                            "type": "string",
                            "enum": ["data_flow", "dependencies", "physics_groups", "category_connections", "all"],
                            "description": "Type of relationships to show (default: 'all')"
                        }
                    }
                }
            )
        ]
    
    @server.call_tool()
    async def call_tool(name: str, arguments: dict):
        if name == "list_actors":
            category = arguments.get("category")
            actors = kb.list_actors(category)
            
            if category:
                response = f"FUSE actors in category '{category}':\n\n"
            else:
                response = "All available FUSE actors:\n\n"
            
            # Group by category if no specific category requested
            if not category:
                for cat in sorted(kb.get_categories()):
                    cat_actors = kb.list_actors(cat)
                    if cat_actors:
                        response += f"**{cat.title()}:**\n"
                        for actor in sorted(cat_actors):
                            response += f"  • {actor}\n"
                        response += "\n"
            else:
                for actor in sorted(actors):
                    actor_info = kb.get_actor(actor)
                    if actor_info:
                        response += f"• **{actor}**: {actor_info.description}\n"
                    else:
                        response += f"• {actor}\n"
            
            return [TextContent(type="text", text=response)]
        
        elif name == "get_actor_info":
            actor_name = arguments["actor_name"]
            actor = kb.get_actor(actor_name)
            
            if not actor:
                return [TextContent(type="text", text=f"❌ Actor '{actor_name}' not found.\n\nTry using 'search_actors' to find the correct name.")]
            
            response = f"# {actor.name} Actor\n\n"
            response += f"**Category**: {actor.category.title()}\n"
            response += f"**Description**: {actor.description}\n\n"
            
            if actor.parameters:
                response += "## Parameters\n\n"
                for param in actor.parameters:
                    response += f"### {param.name}\n"
                    response += f"- **Type**: {param.type}\n"
                    response += f"- **Unit**: {param.unit}\n"
                    response += f"- **Description**: {param.description}\n"
                    if param.default is not None:
                        response += f"- **Default**: {param.default}\n"
                    response += f"- **Required**: {'Yes' if param.required else 'No'}\n\n"
            
            if actor.inputs:
                response += f"**Inputs**: {', '.join(actor.inputs)}\n\n"
            
            if actor.outputs:
                response += f"**Outputs**: {', '.join(actor.outputs)}\n\n"
            
            return [TextContent(type="text", text=response)]
        
        elif name == "search_actors":
            query = arguments["query"]
            results = kb.search_actors(query)
            
            if not results:
                return [TextContent(type="text", text=f"❌ No actors found matching '{query}'.\n\nTry broader terms like 'transport', 'equilibrium', or 'heating'.")]
            
            response = f"🔍 FUSE actors matching '{query}':\n\n"
            for actor in results:
                response += f"**{actor.name}** ({actor.category})\n"
                response += f"  {actor.description}\n\n"
            
            return [TextContent(type="text", text=response)]
        
        elif name == "get_actor_categories":
            categories = kb.get_categories()
            response = "📁 Available FUSE actor categories:\n\n"
            
            for category in sorted(categories):
                actors_in_cat = kb.list_actors(category)
                response += f"• **{category.title()}** ({len(actors_in_cat)} actors)\n"
            
            response += f"\n🎯 Total: {len(kb.actors)} actors across {len(categories)} categories"
            return [TextContent(type="text", text=response)]
        
        elif name == "get_actor_usage_guide":
            actor_name = arguments["actor_name"]
            guide = kb.get_actor_usage_guide(actor_name)
            return [TextContent(type="text", text=guide)]
        
        elif name == "get_actor_relationships":
            return await handle_relationships_tool(arguments)
        
        return [TextContent(type="text", text=f"❌ Unknown tool: {name}")]
    
    @server.list_resources()
    async def list_resources():
        return [
            Resource(
                uri="prompt://fuse_expert",
                name="FUSE Actor Expert System Prompt",
                description="Expert system prompt for FUSE plasma simulation actors with physics guidance",
                mimeType="text/plain"
            ),
            Resource(
                uri="knowledge://fuse_actors",
                name="FUSE Actors Knowledge Base", 
                description="Complete JSON knowledge base of FUSE plasma simulation actors",
                mimeType="application/json"
            )
        ]
    
    @server.read_resource()
    async def read_resource(uri: str):
        if uri == "prompt://fuse_expert":
            prompt_file = Path(__file__).parent / "prompts" / "actor_expert.txt"
            if prompt_file.exists():
                return prompt_file.read_text()
            else:
                return "You are a FUSE plasma simulation expert. Help users understand and work with FUSE actors."
        
        elif uri == "knowledge://fuse_actors":
            import json
            # Create a JSON representation of the knowledge base
            actors_data = {}
            for actor_name, actor_info in kb.actors.items():
                actors_data[actor_name] = {
                    "name": actor_info.name,
                    "category": actor_info.category,
                    "description": actor_info.description,
                    "parameters": [
                        {
                            "name": p.name,
                            "type": p.type,
                            "unit": p.unit,
                            "description": p.description,
                            "default": p.default,
                            "required": p.required
                        } for p in actor_info.parameters
                    ],
                    "inputs": actor_info.inputs,
                    "outputs": actor_info.outputs
                }
            
            return json.dumps({
                "total_actors": len(kb.actors),
                "categories": list(kb.categories.keys()),
                "actors": actors_data
            }, indent=2)
        
        else:
            raise ValueError(f"Unknown resource URI: {uri}")
    
    # Run server
    options = InitializationOptions(
        server_name="fuse-mcp",
        server_version="0.1.0",
        capabilities={"tools": {}, "resources": {}}
    )
    
    async with stdio_server() as (read_stream, write_stream):
        await server.run(read_stream, write_stream, options)


if __name__ == "__main__":
    asyncio.run(main())