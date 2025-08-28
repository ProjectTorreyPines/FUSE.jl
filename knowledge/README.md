# FUSE Actor Knowledge Extractor

Simple Julia tool to extract structured knowledge from FUSE actors using Claude analysis.

## Overview

1. **JSON Schema** (`actor_schema.json`): Defines structure for actor analysis
2. **Simple Extractor** (`extract_knowledge.jl`): Walks actors/ folder, uses Claude to analyze each file
3. **Individual JSONs**: Creates `actors/SOMETHING_actor.json` mirroring source structure  
4. **Combined Knowledge Base**: Single `fuse_knowledge_base.json` with all actors

## Usage

```julia
using Pkg
Pkg.activate(".")

include("extract_knowledge.jl")

# Parallel extraction (recommended - ~6x faster)
kb = extract_fuse_knowledge_parallel()

# Sequential extraction
kb = extract_fuse_knowledge()

# Test with a few actors
test_extraction()
```

## Output Structure

```
knowledge/
├── actor_schema.json           # JSON schema definition
├── fuse_knowledge_base.json    # Combined knowledge base (see below)
└── actors/                     # Individual actor JSONs
    ├── equilibrium/
    │   └── stationary_actor.json
    ├── transport/
    │   └── transport_actor.json
    └── ...
```

### Individual vs Combined Files

**Individual Actor JSONs** (`actors/category/actor_name.json`):
- **Use for**: Focused queries about specific actors
- **Contains**: Detailed actor analysis per schema
- **Benefits**: Easy to update incrementally, mirrors source structure
- **Perfect for**: Development, specific actor documentation

**Combined Knowledge Base** (`fuse_knowledge_base.json`):
- **Use for**: Global queries, MCP server loading, system overview
- **Contains**: All individual actors PLUS:
  - Global metadata and statistics
  - Category-to-actors mappings  
  - Hierarchy classifications (single vs compound)
  - Cross-actor relationships
- **Benefits**: Single file load, fast category queries, system-wide analysis
- **Perfect for**: MCP servers, global search, system architecture overview

## Schema

Each actor JSON contains:
- `name`: Actor class name
- `category`: Domain (equilibrium, transport, etc.) 
- `hierarchy`: "single" or "compound"
- `description`: What the actor does
- `physics_domain`: Physics/engineering domain
- `data_inputs`/`data_outputs`: IMAS data paths
- `sub_actors`: For compound actors (empty for single)
- `key_parameters`: Important parameters
- `usage_notes`: Warnings, best practices

## Requirements

- Julia with JSON package
- Claude Code CLI (`claude` command available)
- FUSE source code at `../src/actors/`

## How It Works

1. Finds all `*_actor.jl` files in `../src/actors/`
2. **Parallel Processing**: Processes actors in batches (8 concurrent Claude requests)
3. For each file: asks Claude to analyze and return JSON per schema
4. Saves individual JSON files mirroring directory structure
5. Combines all into final knowledge base with additional metadata

## Performance

- **Sequential**: ~30+ minutes for 60 actors
- **Parallel (batch_size=8)**: ~3-4 minutes for 60 actors 
- **~8x speed improvement** through batched parallel processing
- Uses Julia threading for concurrent Claude API calls

## Error Handling

The system is designed to be robust against Claude API failures:

1. **Retry Logic**: Each actor gets 1 automatic retry if Claude analysis fails
2. **Individual File Preservation**: Even if Claude fails, any successfully created individual JSONs are preserved
3. **Automatic Recovery**: After extraction, missing actors are automatically recovered from individual files
4. **Graceful Degradation**: System continues processing other actors even if some fail

**Example**: If 2 actors fail Claude analysis but their individual JSONs exist, the system will:
- Complete extraction of remaining 58 actors
- Automatically detect and recover the 2 missing actors from individual files  
- Rebuild relationships with the complete 60-actor dataset
- Report final success: "🔧 Recovered 2 actors from individual files"

Simple, fast, and reliable! 🎯