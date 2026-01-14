#!/bin/bash
# Start FUSE MCP Server with proper environment

# Activate conda environment
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate fuse_mcp

# Change to server directory
cd "$HOME/.julia/dev/FUSE/knowledge/mcp_server/"

# Start the MCP server
exec python server.py
