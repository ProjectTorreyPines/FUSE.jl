#!/bin/zsh
# Start FUSE MCP Server with proper environment
echo "change this file to get it to work:: $(realpath "${BASH_SOURCE[0]}")"

# Activate conda environment
# Use system Python with MCP packages already installed

# Change to server directory
cd "$HOME/.julia/dev/FUSE/knowledge/mcp_server/"

# Start the MCP server
exec python server.py
