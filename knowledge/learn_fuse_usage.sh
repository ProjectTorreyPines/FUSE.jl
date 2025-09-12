#!/bin/bash

# Simple FUSE Learning Script
# Usage: ./learn_fuse_usage.sh notebook1.ipynb notebook2.ipynb ...

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check dependencies
if ! command -v jq &> /dev/null; then
    print_error "jq is required but not found. Please install jq."
    exit 1
fi

if ! command -v claude &> /dev/null; then
    print_error "claude command not found. Please install Claude CLI."
    exit 1
fi

# Filter notebook to extract only input cells
filter_notebook() {
    local notebook="$1"
    local temp_filtered=$(mktemp)
    
    jq '{
        cells: [
            .cells[] | 
            select(.cell_type == "code" or .cell_type == "markdown") |
            {
                cell_type: .cell_type,
                source: .source
            }
        ]
    }' "$notebook" > "$temp_filtered"
    
    echo "$temp_filtered"
}

# Process a single notebook
process_notebook() {
    local notebook="$1"
    local index="$2"
    local total="$3"
    
    print_status "[$index/$total] Processing: $notebook"
    
    if [ ! -f "$notebook" ]; then
        print_error "Notebook not found: $notebook"
        return 1
    fi
    
    # Filter notebook
    print_status "Filtering notebook to extract input cells only..."
    local filtered_notebook
    filtered_notebook=$(filter_notebook "$notebook")
    
    local filtered_size=$(wc -c < "$filtered_notebook")
    local original_size=$(wc -c < "$notebook")
    local reduction_pct=$(( 100 - (filtered_size * 100 / original_size) ))
    print_status "Size reduction: ${original_size} → ${filtered_size} bytes (${reduction_pct}% smaller)"
    
    # Create prompt
    local prompt="I am learning FUSE by studying Jupyter notebooks.

**Current Notebook:** $notebook (filtered to show only input cells)

**Task:**
Study this filtered notebook and update how_to_use_fuse.md with new FUSE knowledge you learn.

The notebook contains only code and markdown input cells - no outputs or execution results.

Focus on:
- FUSE usage patterns and workflows
- Code examples that show how to use FUSE
- Important concepts and best practices
- Parameter configurations and their effects

Build incrementally on existing knowledge in how_to_use_fuse.md

Keep a list of tutorial files at the end of how_to_use_fuse.md that were used to learn,
and use this information to be aware that the user may have already given you the same notebook in the past.

**Filtered notebook to process:** $filtered_notebook"
    
    # Process with Claude
    print_status "Processing with Claude..."
    if claude --allowed-tools "Read Edit Write" --print <<< "$prompt"; then
        print_success "Processed: $notebook"
        rm -f "$filtered_notebook"
        return 0
    else
        print_error "Failed to process: $notebook"
        rm -f "$filtered_notebook"
        return 1
    fi
}

# Main function
main() {
    if [ $# -eq 0 ]; then
        echo "Usage: $0 notebook1.ipynb notebook2.ipynb ..."
        echo ""
        echo "Example:"
        echo "  $0 examples/tutorial.ipynb examples/tutorial_imas.ipynb"
        echo ""
        echo "This script will:"
        echo "  1. Filter each notebook to extract only input cells"
        echo "  2. Process each notebook with Claude"
        echo "  3. Update how_to_use_fuse.md with learned knowledge"
        exit 1
    fi
    
    print_status "Starting simple FUSE learning process..."
    print_status "Notebooks to process: $#"
    
    local processed=0
    local failed=0
    local index=0
    
    for notebook in "$@"; do
        index=$((index + 1))
        if process_notebook "$notebook" "$index" "$#"; then
            processed=$((processed + 1))
        else
            failed=$((failed + 1))
        fi
        
        # Small delay between notebooks
        sleep 2
    done
    
    print_success "Learning completed!"
    print_status "Processed: $processed notebooks"
    if [ $failed -gt 0 ]; then
        print_error "Failed: $failed notebooks"
    fi
    print_status "Check how_to_use_fuse.md for updated knowledge"
}

# Run main function
main "$@"