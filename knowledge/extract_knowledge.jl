"""
FUSE Actor Knowledge Extractor

Simple approach:
1. Walk actors/ folder to find *_actor.jl files
2. Use Claude to analyze each file and return JSON according to schema
3. Save individual JSON files mirroring directory structure
4. Generate combined knowledge base
"""

using JSON
using Dates

"""
    extract_all_actors(src_dir::String, knowledge_dir::String="knowledge")

Extract knowledge for all FUSE actors using Claude analysis.
Creates individual JSON files and combined knowledge base.
"""
function extract_all_actors(src_dir::String, knowledge_dir::String="knowledge"; limit=nothing)
    actors_src_dir = joinpath(src_dir, "actors")
    actors_kb_dir = joinpath(knowledge_dir, "actors")
    
    if !isdir(actors_src_dir)
        error("Actors source directory not found: $actors_src_dir")
    end
    
    println("🔍 Extracting FUSE Actor Knowledge...")
    println("Source: $actors_src_dir")
    println("Output: $actors_kb_dir")
    
    # Ensure output directory exists
    mkpath(actors_kb_dir)
    
    # Load schema for validation
    schema_path = joinpath(knowledge_dir, "actor_schema.json")
    if !isfile(schema_path)
        error("Schema file not found: $schema_path")
    end
    
    # Find all actor files
    actor_files = find_actor_files(actors_src_dir)
    if limit !== nothing
        actor_files = actor_files[1:min(limit, length(actor_files))]
    end
    println("Found $(length(actor_files)) actor files")
    
    all_actors = Dict()
    categories = Dict()
    
    # Process each actor
    for (i, (actor_file, category)) in enumerate(actor_files)
        println("[$i/$(length(actor_files))] Processing $(basename(actor_file))...")
        
        try
            # Analyze with Claude
            actor_data = analyze_actor_with_claude(actor_file, category, schema_path)
            
            if haskey(actor_data, "error")
                @warn "Failed to analyze $(basename(actor_file)): $(actor_data["error"])"
                continue
            end
            
            # Save individual JSON file
            save_individual_actor_json(actor_data, actor_file, actors_src_dir, actors_kb_dir)
            
            # Track for combined knowledge base
            actor_name = actor_data["name"]
            all_actors[actor_name] = actor_data
            
            # Track categories
            if !haskey(categories, category)
                categories[category] = []
            end
            push!(categories[category], actor_name)
            
        catch e
            @warn "Error processing $(basename(actor_file)): $e"
        end
        
        # Small delay to be nice to Claude
        sleep(0.5)
    end
    
    # Generate combined knowledge base
    knowledge_base = create_knowledge_base(all_actors, categories)
    
    # Save combined knowledge base
    kb_path = joinpath(knowledge_dir, "fuse_knowledge_base.json")
    open(kb_path, "w") do io
        JSON.print(io, knowledge_base, 2)
    end
    
    println("✅ Extraction complete!")
    println("Individual actors: $actors_kb_dir")
    println("Combined knowledge base: $kb_path")
    println("Processed $(length(all_actors)) actors across $(length(categories)) categories")
    
    return knowledge_base
end

"""
    find_actor_files(actors_dir::String) -> Vector{Tuple{String, String}}

Find all *_actor.jl files and determine their categories from directory structure.
Returns vector of (filepath, category) tuples.
"""
function find_actor_files(actors_dir::String)
    actor_files = []
    
    for (root, dirs, files) in walkdir(actors_dir)
        # Determine category from directory structure
        category = basename(root)
        if category == "actors"
            category = "uncategorized"  # Root actors directory
        end
        
        for file in files
            if endswith(file, "_actor.jl")
                filepath = joinpath(root, file)
                push!(actor_files, (filepath, category))
            end
        end
    end
    
    return actor_files
end

"""
    analyze_actor_with_claude(actor_file::String, category::String, schema_path::String) -> Dict

Use Claude to analyze a single actor file and return structured JSON.
"""
function analyze_actor_with_claude(actor_file::String, category::String, schema_path::String)
    # Read actor source code
    if !isfile(actor_file)
        return Dict("error" => "File not found: $actor_file")
    end
    
    source_code = read(actor_file, String)
    
    # Read schema for Claude
    schema_content = read(schema_path, String)
    
    # Extract actor name from filename
    filename = basename(actor_file)
    actor_name = "Actor" * titlecase(replace(filename, "_actor.jl" => ""))
    
    # Create analysis prompt
    prompt = """
You are analyzing a FUSE (fusion simulation) actor. Please analyze this Julia source code and return a JSON object that exactly matches the provided schema.

**CRITICAL INSTRUCTIONS:**
- Return ONLY valid JSON, no markdown, no explanations, no code blocks
- The JSON must validate against the schema
- Fill out ALL required fields
- Use the category "$category" for the category field
- Extract the actual actor name from the struct definition in the code
- For sub_actors: only include if this is a compound actor (inherits from CompoundAbstractActor)
- For data_inputs/data_outputs: look for IMAS data paths in the code
- Be concise but informative in descriptions

**JSON Schema:**
```json
$schema_content
```

**Actor Source Code:**
```julia
$source_code
```

Return the JSON object:"""

    try
        # Use Claude Code to analyze
        result = read(`claude "$prompt"`, String)
        
        # Clean up result - remove any markdown formatting
        json_str = strip(result)
        
        # Remove markdown code blocks if present
        if startswith(json_str, "```json")
            json_str = replace(json_str, r"```json\s*" => "")
            json_str = replace(json_str, r"\s*```\s*$" => "")
        elseif startswith(json_str, "```")
            json_str = replace(json_str, r"```\s*" => "")
            json_str = replace(json_str, r"\s*```\s*$" => "")
        end
        
        # Parse JSON
        actor_data = JSON.parse(json_str)
        
        # Validate required fields exist
        required_fields = ["name", "category", "hierarchy", "description", "physics_domain", 
                          "data_inputs", "data_outputs", "sub_actors", "key_parameters", "usage_notes"]
        for field in required_fields
            if !haskey(actor_data, field)
                return Dict("error" => "Missing required field: $field")
            end
        end
        
        return actor_data
        
    catch e
        return Dict("error" => "Claude analysis failed: $e")
    end
end

"""
    save_individual_actor_json(actor_data::Dict, original_file::String, src_base::String, output_base::String)

Save individual actor JSON file maintaining directory structure.
"""
function save_individual_actor_json(actor_data::Dict, original_file::String, src_base::String, output_base::String)
    # Calculate relative path from source
    rel_path = relpath(original_file, src_base)
    
    # Replace .jl with .json
    json_filename = replace(basename(rel_path), ".jl" => ".json")
    json_rel_path = joinpath(dirname(rel_path), json_filename)
    
    # Full output path
    output_path = joinpath(output_base, json_rel_path)
    
    # Ensure directory exists
    mkpath(dirname(output_path))
    
    # Save JSON
    open(output_path, "w") do io
        JSON.print(io, actor_data, 2)
    end
end

"""
    create_knowledge_base(all_actors::Dict, categories::Dict) -> Dict

Create combined knowledge base from individual actor analyses.
"""
function create_knowledge_base(all_actors::Dict, categories::Dict)
    single_count = count(actor -> actor["hierarchy"] == "single", values(all_actors))
    compound_count = count(actor -> actor["hierarchy"] == "compound", values(all_actors))
    
    return Dict(
        "metadata" => Dict(
            "generated_at" => string(now()),
            "extraction_method" => "claude_analysis",
            "total_actors" => length(all_actors),
            "single_actors" => single_count,
            "compound_actors" => compound_count,
            "categories" => length(categories),
            "version" => "1.0.0"
        ),
        "actors" => all_actors,
        "categories" => categories,
        "hierarchy" => Dict(
            "single" => [name for (name, actor) in all_actors if actor["hierarchy"] == "single"],
            "compound" => [name for (name, actor) in all_actors if actor["hierarchy"] == "compound"]
        )
    )
end

# Convenience function
"""
    extract_fuse_knowledge()

Extract knowledge using default paths (assumes running from knowledge/ directory).
"""
# Test function for a few actors
function test_extraction()
    extract_all_actors("../src", ".", limit=3)
end

function extract_fuse_knowledge()
    extract_all_actors("../src", ".")
end