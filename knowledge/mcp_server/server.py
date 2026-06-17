#!/usr/bin/env python3
"""
FUSE MCP Server - Model Context Protocol server for FUSE plasma simulation actors

This server provides LLM access to FUSE actor documentation, usage guides, parameter information,
relationships between actors, and use case definitions. It supports comprehensive queries about 
plasma simulation workflows and tokamak configurations.
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
usage_guide_content = None
usage_guide_sections = None

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
            print(f"Loaded relationships: {list(relationships_data.keys())}", file=sys.stderr)
    except Exception as e:
        print(f"Failed to load relationships: {e}", file=sys.stderr)
        relationships_data = {}

def load_usage_guide():
    """Load FUSE usage guide from how_to_use_fuse.md"""
    global usage_guide_content, usage_guide_sections
    import re
    
    server_dir = Path(__file__).parent
    usage_guide_file = server_dir.parent / "how_to_use_fuse.md"
    
    try:
        with open(usage_guide_file, 'r', encoding='utf-8') as f:
            usage_guide_content = f.read()
        
        # Parse sections based on headers
        usage_guide_sections = parse_markdown_sections(usage_guide_content)
        print(f"Loaded usage guide with {len(usage_guide_sections)} sections", file=sys.stderr)

    except Exception as e:
        print(f"Failed to load usage guide: {e}", file=sys.stderr)
        usage_guide_content = ""
        usage_guide_sections = {}

def parse_markdown_sections(content):
    """Parse markdown content into sections based on headers"""
    import re
    
    sections = {}
    
    # Split by ## headers (main sections)
    parts = re.split(r'^## (.+)$', content, flags=re.MULTILINE)
    
    if len(parts) > 1:
        # First part is usually intro/overview
        if parts[0].strip():
            sections["overview"] = parts[0].strip()
        
        # Parse main sections
        for i in range(1, len(parts), 2):
            if i + 1 < len(parts):
                section_title = parts[i].strip()
                section_content = parts[i + 1].strip()
                
                # Create key from title (lowercase, replace spaces with underscores)
                section_key = normalize_topic_key(section_title)
                sections[section_key] = f"## {section_title}\n\n{section_content}"
    
    # Create dynamic aliases based on content analysis
    sections.update(create_dynamic_aliases(sections))
    
    return sections

def normalize_topic_key(title):
    """Normalize section title to create a clean topic key"""
    import re
    return re.sub(r'[^\w\s]', '', title.lower()).replace(' ', '_').replace('__', '_').strip('_')

def create_dynamic_aliases(sections):
    """Create simple aliases from section titles themselves"""
    import re
    aliases = {}
    
    for section_key, section_content in sections.items():
        # Extract title from content
        title_match = re.search(r'^## (.+)$', section_content, re.MULTILINE)
        if title_match:
            title = title_match.group(1)
            
            # Create aliases from individual words in the title
            words = re.findall(r'\b\w+\b', title.lower())
            for word in words:
                if len(word) > 3 and word not in aliases and word != section_key:
                    aliases[word] = section_content
    
    return aliases

def get_usage_guide_topics():
    """Get available topics in the usage guide dynamically from headers"""
    if not usage_guide_sections:
        return []
    
    # Get main sections by finding keys that start with "##" in their content (original sections)
    main_topics = []
    for key, content in usage_guide_sections.items():
        if content.startswith("##") or content.startswith("# "):
            main_topics.append(key)
    
    return sorted(main_topics)

def get_topic_description(topic_key):
    """Extract a brief description from the topic's content"""
    if topic_key not in usage_guide_sections:
        return "Topic information"
    
    content = usage_guide_sections[topic_key]
    
    # Extract first line after header as description
    lines = content.split('\n')
    for i, line in enumerate(lines):
        if line.startswith('##'):
            # Look for first non-empty line after header
            for j in range(i + 1, min(i + 5, len(lines))):
                if lines[j].strip() and not lines[j].startswith('#'):
                    # Take first sentence or first line, whichever is shorter
                    desc = lines[j].strip()
                    if '.' in desc:
                        desc = desc.split('.')[0] + '.'
                    return desc[:100] + ('...' if len(desc) > 100 else '')
    
    # Fallback: use topic key as description
    return f"Information about {topic_key.replace('_', ' ')}"

async def handle_usage_guide_tool(arguments: dict):
    """Handle the get_fuse_usage_guide tool"""
    if usage_guide_sections is None:
        return [TextContent(type="text", text="❌ Usage guide not loaded")]
    
    topic = arguments.get("topic", "").lower()
    
    # If no topic specified, show available topics
    if not topic:
        response = "# 📚 FUSE Usage Guide Topics\n\n"
        response += "Available topics for detailed guidance:\n\n"
        
        # Get main topics dynamically
        main_topics = get_usage_guide_topics()
        
        for topic_key in main_topics:
            description = get_topic_description(topic_key)
            response += f"• **{topic_key}**: {description}\n"
        
        # Show available aliases/search terms
        alias_topics = [key for key in sorted(usage_guide_sections.keys()) if key not in main_topics]
        if alias_topics:
            response += f"\n🔍 **Search terms**: {', '.join(alias_topics[:15])}" # Show first 15 aliases
            if len(alias_topics) > 15:
                response += f", and {len(alias_topics) - 15} more..."
        
        response += "\n\n💡 **Usage**: Use `get_fuse_usage_guide` with a specific topic to get detailed information."
        
        return [TextContent(type="text", text=response)]
    
    # Search for topic
    topic_content = None
    
    # Direct match
    if topic in usage_guide_sections:
        topic_content = usage_guide_sections[topic]
    else:
        # Fuzzy search - find sections containing the topic
        matching_sections = []
        for key, content in usage_guide_sections.items():
            if topic in key or topic in content.lower():
                matching_sections.append((key, content))
        
        if len(matching_sections) == 1:
            topic_content = matching_sections[0][1]
        elif len(matching_sections) > 1:
            response = f"# 🔍 Multiple sections found for '{topic}'\n\n"
            for key, _ in matching_sections:
                response += f"• {key}\n"
            response += f"\nPlease specify one of these topics for detailed information."
            return [TextContent(type="text", text=response)]
    
    if topic_content:
        return [TextContent(type="text", text=topic_content)]
    else:
        available_topics = ", ".join(sorted(usage_guide_sections.keys())[:10])
        return [TextContent(type="text", text=f"❌ Topic '{topic}' not found.\n\nAvailable topics: {available_topics}...\n\nUse `get_fuse_usage_guide` without arguments to see all available topics.")]

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

    print(f"Loading actors from: {actors_dir}", file=sys.stderr)
    print(f"Knowledge base file: {knowledge_file}", file=sys.stderr)
    
    kb = ActorKnowledgeBase(
        actors_dir=actors_dir,
        knowledge_base_file=knowledge_file
    )
    
    # Load relationships data
    load_relationships()
    
    # Load usage guide
    load_usage_guide()
    
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
            ),
            Tool(
                name="get_fuse_usage_guide",
                description="Get FUSE usage guide for workflows, concepts, and best practices. Covers core concepts, transport modeling, time-dependent simulations, IMAS data manipulation, and more",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "topic": {
                            "type": "string",
                            "description": "Specific topic to get guidance on (e.g., 'transport', 'flux_matching', 'time_dependent', 'actors', 'workflows', 'imas'). Leave empty to see all available topics."
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
        
        elif name == "get_fuse_usage_guide":
            return await handle_usage_guide_tool(arguments)
        
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
            ),
            Resource(
                uri="knowledge://fuse_usage_guide",
                name="FUSE Usage Guide",
                description="Comprehensive guide for FUSE workflows, concepts, transport modeling, time-dependent simulations, and best practices",
                mimeType="text/markdown"
            )
        ]
    
    @server.read_resource()
    async def read_resource(uri: str):
        from mcp.types import TextResourceContents

        if uri == "prompt://fuse_expert":
            prompt_file = Path(__file__).parent / "prompts" / "actor_expert.txt"
            if prompt_file.exists():
                content = prompt_file.read_text()
            else:
                content = "You are a FUSE plasma simulation expert. Help users understand and work with FUSE actors."

            return TextResourceContents(
                uri=uri,
                text=content,
                mimeType="text/plain"
            )

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

            content = json.dumps({
                "total_actors": len(kb.actors),
                "categories": list(kb.categories.keys()),
                "actors": actors_data
            }, indent=2)

            return TextResourceContents(
                uri=uri,
                text=content,
                mimeType="application/json"
            )

        elif uri == "knowledge://fuse_usage_guide":
            if usage_guide_content:
                content = usage_guide_content
            else:
                content = "# FUSE Usage Guide\n\nUsage guide not loaded. Please check server configuration."

            return TextResourceContents(
                uri=uri,
                text=content,
                mimeType="text/markdown"
            )

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