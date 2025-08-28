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

# Extract all actors (from knowledge/ directory)
kb = extract_fuse_knowledge()

# Or specify custom paths
kb = extract_all_actors("../src", ".")
```

## Output Structure

```
knowledge/
├── actor_schema.json           # JSON schema definition
├── fuse_knowledge_base.json    # Combined knowledge base
└── actors/                     # Individual actor JSONs
    ├── equilibrium/
    │   └── stationary_actor.json
    ├── transport/
    │   └── transport_actor.json
    └── ...
```

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
2. For each file: asks Claude to analyze and return JSON per schema
3. Saves individual JSON files mirroring directory structure
4. Combines all into final knowledge base

Simple and effective!