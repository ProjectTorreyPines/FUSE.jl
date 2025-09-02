#!/bin/bash
# Start FUSE MCP Server with proper environment

# Activate conda environment
source ~/miniconda3/bin/activate
conda activate fuse_mcp

# Change to server directory
cd "/Users/meneghini/.julia/dev/FUSE/knowledge/mcp_server"

# Start the MCP server
exec python server.py

