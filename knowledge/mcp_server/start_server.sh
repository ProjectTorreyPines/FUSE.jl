#!/bin/zsh
# Start FUSE MCP Server with proper environment
echo "change this file to get it to work:: $(realpath "${BASH_SOURCE[0]}")"

# Activate conda environment
source $HOME/.zshrc
mamba activate fuse_mcp # this is your conda env

# Change to server directory
cd "$HOME/.julia/dev/FUSE/knowledge/mcp_server/"

# Start the MCP server
exec python server.py
