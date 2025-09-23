#!/bin/zsh
# Start FUSE MCP Server with proper environment
echo "change this file to get it to work:: $(realpath "${BASH_SOURCE[0]}")"

# Activate conda environment
<<<<<<< HEAD
source $HOME/.zshrc
conda activate mcp_server # this is your conda env

# Change to server directory
cd "$HOME/.julia/dev/FUSE/knowledge/mcp_server/"
=======
source /Users/tims/.zshrc
conda activate mcp_server # this is your conda env

# Change to server directory
cd "/Users/tims/.julia/dev/FUSE/knowledge/mcp_server/"
>>>>>>> parent of 4865739a (Revert "Merge branch 'master' into fix/optional-pedestal-density-tanh")

# Start the MCP server
exec python server.py
