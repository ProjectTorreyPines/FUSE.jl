# FUSE MCP Server

Model Context Protocol (MCP) server providing LLM access to FUSE plasma simulation actors.

## Features

- **Actor Knowledge Base**: Access to complete FUSE actor documentation
- **Parameter Guidance**: Detailed parameter descriptions with units and physics context  
- **Usage Examples**: Actor usage patterns and best practices
- **Search & Discovery**: Find actors by category, name, or physics domain
- **Expert System**: Physics-aware guidance for plasma simulation workflows

## Installation

```bash
cd knowledge/mcp_server
pip install -r requirements.txt
```

## Usage

### Running the Server

```bash
./start_server.sh
# or directly:
python server.py
```

### Available Tools

- `list_actors` - List all available actors, optionally by category
- `get_actor_info` - Get detailed information about a specific actor
- `search_actors` - Search actors by name, description, or physics domain
- `get_actor_categories` - List all actor categories
- `get_actor_usage_guide` - Get comprehensive usage guide for an actor

### Available Resources

- `prompt://actor_expert` - Expert system prompt for FUSE physics guidance
- `knowledge://fuse_actors` - Complete actors knowledge base (JSON)

## Architecture

```
mcp_server/
├── server.py               # Main MCP server
├── start_server.sh         # Startup script
├── requirements.txt        # Python dependencies
├── fuse_mcp/              # Actor knowledge module
│   ├── __init__.py
│   ├── actor_knowledge.py # Knowledge base management
│   └── tools/
│       ├── __init__.py
│       └── actor_info.py  # Actor information tools
└── prompts/
    └── actor_expert.txt    # Physics expert prompt
```

## Knowledge Base

The server loads actor information from:
- Individual JSON files in `../actors/` (organized by category)
- Unified knowledge base `../fuse_knowledge_base.json`
- Actor parameters, inputs, outputs, and usage examples

## Future Enhancements (Phase 2)

- **Actor Execution**: Direct actor execution via Julia RemoteREPL
- **Result Analysis**: Parse and interpret actor outputs  
- **Workflow Orchestration**: Multi-actor simulation workflows
- **Parameter Validation**: Check parameter compatibility and physics constraints

## Integration

Configure in Claude Code or other MCP clients:

```json
{
  "mcpServers": {
    "fuse": {
      "command": "python",
      "args": ["/path/to/fuse/knowledge/mcp_server/server.py"]
    }
  }
}
```