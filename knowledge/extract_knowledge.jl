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
using Base.Threads

"""
    extract_all_actors(src_dir::String, knowledge_dir::String="knowledge")

Extract knowledge for all FUSE actors using Claude analysis.
Creates individual JSON files and combined knowledge base.
"""
function extract_all_actors_parallel(src_dir::String, knowledge_dir::String="knowledge"; batch_size=8, limit=nothing)
    actors_src_dir = joinpath(src_dir, "actors")
    actors_kb_dir = joinpath(knowledge_dir, "actors")
    
    if !isdir(actors_src_dir)
        error("Actors source directory not found: $actors_src_dir")
    end
    
    println("🚀 Extracting FUSE Actor Knowledge (Parallel Batches)...")
    println("Source: $actors_src_dir")
    println("Output: $actors_kb_dir")
    println("Batch size: $batch_size concurrent requests")
    println("Available threads: $(Threads.nthreads())")
    
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
    
    # Create batches
    batches = [actor_files[i:min(i+batch_size-1, end)] 
              for i in 1:batch_size:length(actor_files)]
    println("Processing in $(length(batches)) batches")
    
    # Process batches
    total_processed = 0
    for (batch_num, batch) in enumerate(batches)
        batch_start = total_processed + 1
        batch_end = total_processed + length(batch)
        println("\n📦 Batch $batch_num/$(length(batches)): Processing actors [$batch_start-$batch_end]")
        
        # Start parallel tasks for this batch
        tasks = []
        for (actor_file, category) in batch
            task = Threads.@spawn begin
                try
                    analyze_actor_with_claude_retry(actor_file, category, schema_path)
                catch e
                    Dict("error" => "Task failed: $e", "file" => actor_file)
                end
            end
            push!(tasks, (task, actor_file, category))
        end
        
        # Wait for all tasks in batch to complete and collect results
        batch_results = []
        for (i, (task, actor_file, category)) in enumerate(tasks)
            actor_data = fetch(task)
            push!(batch_results, (actor_data, actor_file, category))
            print("✓")  # Progress indicator
            flush(stdout)
        end
        
        # Process results from this batch
        for (actor_data, actor_file, category) in batch_results
            if haskey(actor_data, "error")
                @warn "Failed to analyze $(basename(actor_file)): $(actor_data["error"])"
                continue
            end
            
            try
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
                @warn "Error saving $(basename(actor_file)): $e"
            end
        end
        
        total_processed += length(batch)
        println("\n   Completed: $total_processed/$(length(actor_files)) actors")
        
        # Small delay between batches to be nice to Claude API
        if batch_num < length(batches)
            sleep(1.0)
        end
    end
    
    # Analyze relationships between actors
    relationships = analyze_actor_relationships(all_actors)
    
    # Generate combined knowledge base
    knowledge_base = create_knowledge_base(all_actors, categories, relationships)
    
    # Save combined knowledge base
    kb_path = joinpath(knowledge_dir, "fuse_knowledge_base.json")
    open(kb_path, "w") do io
        JSON.print(io, knowledge_base, 2)
    end
    
    # Check for and recover any missing actors from individual files
    missing_actors = recover_missing_actors!(all_actors, categories, actors_kb_dir)
    
    # Check for orphaned individual files (JSONs without corresponding source files)
    orphaned_actors = check_orphaned_actors(actors_kb_dir, src_dir)
    
    if !isempty(missing_actors) || !isempty(orphaned_actors)
        if !isempty(missing_actors)
            println("🔧 Recovered $(length(missing_actors)) actors from individual files: $(join(missing_actors, ", "))")
        end
        if !isempty(orphaned_actors)
            println("⚠️  Found $(length(orphaned_actors)) orphaned actor files (source deleted/moved): $(join(orphaned_actors, ", "))")
        end
        
        # Re-analyze relationships with complete actor set
        relationships = analyze_actor_relationships(all_actors)
        knowledge_base = create_knowledge_base(all_actors, categories, relationships)
        
        # Re-save with complete data
        open(kb_path, "w") do io
            JSON.print(io, knowledge_base, 2)
        end
    end
    
    println("\n✅ Parallel extraction complete!")
    println("Individual actors: $actors_kb_dir")
    println("Combined knowledge base: $kb_path")
    println("Processed $(length(all_actors)) actors across $(length(categories)) categories")
    
    return knowledge_base
end

# Keep original sequential version for comparison
function extract_all_actors(src_dir::String, knowledge_dir::String="knowledge"; limit=nothing)
    return extract_all_actors_parallel(src_dir, knowledge_dir; batch_size=1, limit=limit)
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
    analyze_actor_with_claude_retry(actor_file::String, category::String, schema_path::String) -> Dict

Use Claude to analyze a single actor file with retry logic.
"""
function analyze_actor_with_claude_retry(actor_file::String, category::String, schema_path::String; max_retries=1)
    for attempt in 1:(max_retries + 1)
        result = analyze_actor_with_claude(actor_file, category, schema_path)
        
        if !haskey(result, "error")
            return result
        end
        
        if attempt <= max_retries
            @info "Retrying Claude analysis for $(basename(actor_file)) (attempt $(attempt + 1))"
            sleep(1.0)  # Brief pause before retry
        end
    end
    
    # All retries failed
    return Dict("error" => "Claude analysis failed after $(max_retries + 1) attempts")
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
    analyze_actor_relationships(all_actors::Dict) -> Dict

Analyze relationships between actors based on their data flows, hierarchies, and domains.
"""
function analyze_actor_relationships(all_actors::Dict)
    println("🔗 Analyzing actor relationships...")
    
    relationships = Dict(
        "data_flow_connections" => analyze_data_flow_connections(all_actors),
        "sub_actor_dependencies" => analyze_sub_actor_dependencies(all_actors), 
        "physics_domain_groups" => analyze_physics_domain_groups(all_actors),
        "category_connections" => analyze_category_connections(all_actors)
    )
    
    return relationships
end

"""
    analyze_data_flow_connections(all_actors::Dict) -> Dict

Find data flow connections where Actor A's outputs match Actor B's inputs.
"""
function analyze_data_flow_connections(all_actors::Dict)
    connections = Dict()
    
    for (actor_a, info_a) in all_actors
        connections[actor_a] = Dict("feeds_to" => [], "receives_from" => [])
        
        for (actor_b, info_b) in all_actors
            if actor_a == actor_b
                continue
            end
            
            # Check if A's outputs match B's inputs
            outputs_a = Set(info_a["data_outputs"])
            inputs_b = Set(info_b["data_inputs"])
            
            # Find overlapping data paths (exact matches or prefix matches)
            shared_data = []
            for output in outputs_a
                for input in inputs_b
                    if output == input || startswith(input, output) || startswith(output, input)
                        push!(shared_data, Dict("data_path" => output, "connection_type" => "data_flow"))
                    end
                end
            end
            
            if !isempty(shared_data)
                push!(connections[actor_a]["feeds_to"], Dict("actor" => actor_b, "shared_data" => shared_data))
                if !haskey(connections, actor_b)
                    connections[actor_b] = Dict("feeds_to" => [], "receives_from" => [])
                end
                push!(connections[actor_b]["receives_from"], Dict("actor" => actor_a, "shared_data" => shared_data))
            end
        end
    end
    
    return connections
end

"""
    analyze_sub_actor_dependencies(all_actors::Dict) -> Dict

Analyze compound actor → sub-actor relationships.
"""
function analyze_sub_actor_dependencies(all_actors::Dict)
    dependencies = Dict("compound_to_sub" => Dict(), "sub_to_compound" => Dict())
    
    for (actor_name, info) in all_actors
        if info["hierarchy"] == "compound" && !isempty(info["sub_actors"])
            dependencies["compound_to_sub"][actor_name] = [sa["name"] for sa in info["sub_actors"]]
            
            # Build reverse mapping - only for sub-actors that exist in our extracted set
            for sub_actor in info["sub_actors"]
                sub_name = sub_actor["name"]
                
                # Check if this sub-actor actually exists in our extracted actors
                if haskey(all_actors, sub_name)
                    if !haskey(dependencies["sub_to_compound"], sub_name)
                        dependencies["sub_to_compound"][sub_name] = []
                    end
                    push!(dependencies["sub_to_compound"][sub_name], actor_name)
                else
                    # Note: Some sub-actors might not be extracted (they might not follow *_actor.jl naming)
                    @debug "Sub-actor $sub_name referenced by $actor_name not found in extracted actors"
                end
            end
        end
    end
    
    return dependencies
end

"""
    analyze_physics_domain_groups(all_actors::Dict) -> Dict

Group actors by physics domain for domain-based relationships.
"""
function analyze_physics_domain_groups(all_actors::Dict)
    domain_groups = Dict()
    
    for (actor_name, info) in all_actors
        domain = info["physics_domain"]
        if !haskey(domain_groups, domain)
            domain_groups[domain] = []
        end
        push!(domain_groups[domain], actor_name)
    end
    
    return domain_groups
end

"""
    analyze_category_connections(all_actors::Dict) -> Dict

Analyze connections between different categories of actors.
"""
function analyze_category_connections(all_actors::Dict)
    category_connections = Dict()
    
    # Group actors by category
    by_category = Dict()
    for (actor_name, info) in all_actors
        category = info["category"]
        if !haskey(by_category, category)
            by_category[category] = []
        end
        push!(by_category[category], actor_name)
    end
    
    # Find inter-category data flow connections
    for (cat_a, actors_a) in by_category
        category_connections[cat_a] = Dict()
        
        for (cat_b, actors_b) in by_category
            if cat_a == cat_b
                continue
            end
            
            connections = 0
            for actor_a in actors_a
                for actor_b in actors_b
                    info_a = all_actors[actor_a]
                    info_b = all_actors[actor_b]
                    
                    # Check for data flow connections
                    outputs_a = Set(info_a["data_outputs"])
                    inputs_b = Set(info_b["data_inputs"])
                    
                    if !isempty(intersect(outputs_a, inputs_b))
                        connections += 1
                    end
                end
            end
            
            if connections > 0
                category_connections[cat_a][cat_b] = connections
            end
        end
    end
    
    return category_connections
end

"""
    recover_missing_actors!(all_actors::Dict, categories::Dict, actors_kb_dir::String) -> Vector{String}

Check for actors that exist as individual JSON files but are missing from the combined knowledge base.
Add any missing actors and return their names.
"""
function recover_missing_actors!(all_actors::Dict, categories::Dict, actors_kb_dir::String)
    missing_actors = []
    
    # Walk through all individual JSON files
    for (root, dirs, files) in walkdir(actors_kb_dir)
        category = basename(root)
        if category == "actors"
            category = "uncategorized"
        end
        
        for file in files
            if endswith(file, ".json")
                filepath = joinpath(root, file)
                try
                    actor_data = JSON.parsefile(filepath)
                    actor_name = actor_data["name"]
                    
                    # Check if this actor is missing from the combined knowledge base
                    if !haskey(all_actors, actor_name)
                        @info "Recovering missing actor: $actor_name from $filepath"
                        
                        # Add to actors
                        all_actors[actor_name] = actor_data
                        
                        # Add to categories
                        actor_category = actor_data["category"]
                        if !haskey(categories, actor_category)
                            categories[actor_category] = []
                        end
                        push!(categories[actor_category], actor_name)
                        
                        push!(missing_actors, actor_name)
                    end
                catch e
                    @warn "Failed to parse individual JSON file $filepath: $e"
                end
            end
        end
    end
    
    return missing_actors
end

"""
    check_orphaned_actors(actors_kb_dir::String, src_dir::String) -> Vector{String}

Find individual JSON files that no longer have corresponding source files.
This handles cases where source files were moved, renamed, or deleted.
"""
function check_orphaned_actors(actors_kb_dir::String, src_dir::String)
    orphaned_actors = []
    actors_src_dir = joinpath(src_dir, "actors")
    
    # Get all existing source file names (without path)
    existing_source_files = Set{String}()
    if isdir(actors_src_dir)
        for (root, dirs, files) in walkdir(actors_src_dir)
            for file in files
                if endswith(file, "_actor.jl")
                    push!(existing_source_files, file)
                end
            end
        end
    end
    
    # Check each individual JSON file
    for (root, dirs, files) in walkdir(actors_kb_dir)
        for file in files
            if endswith(file, ".json")
                # Convert JSON filename to expected source filename
                expected_source = replace(file, ".json" => ".jl")
                
                if !in(expected_source, existing_source_files)
                    filepath = joinpath(root, file)
                    try
                        actor_data = JSON.parsefile(filepath)
                        actor_name = actor_data["name"]
                        push!(orphaned_actors, actor_name)
                        @warn "Orphaned actor JSON found: $actor_name (source file $expected_source not found)"
                    catch e
                        @warn "Could not parse potentially orphaned JSON file $filepath: $e"
                    end
                end
            end
        end
    end
    
    return orphaned_actors
end

"""
    clean_orphaned_actors(actors_kb_dir::String, src_dir::String; interactive=true) -> Int

Remove orphaned individual JSON files that no longer have corresponding source files.
Returns the number of files cleaned.
"""
function clean_orphaned_actors(actors_kb_dir::String, src_dir::String; interactive=true)
    orphaned_actors = check_orphaned_actors(actors_kb_dir, src_dir)
    
    if isempty(orphaned_actors)
        println("No orphaned actor files found.")
        return 0
    end
    
    println("Found $(length(orphaned_actors)) orphaned actor files:")
    for actor in orphaned_actors
        println("  - $actor")
    end
    
    if interactive
        println("\nRemove these orphaned files? (y/N): ")
        response = strip(readline())
        if lowercase(response) != "y"
            println("Cleanup cancelled.")
            return 0
        end
    end
    
    # Remove orphaned files
    removed_count = 0
    for (root, dirs, files) in walkdir(actors_kb_dir)
        for file in files
            if endswith(file, ".json")
                filepath = joinpath(root, file)
                try
                    actor_data = JSON.parsefile(filepath)
                    actor_name = actor_data["name"]
                    
                    if actor_name in orphaned_actors
                        rm(filepath)
                        println("🗑️  Removed orphaned file: $filepath")
                        removed_count += 1
                    end
                catch e
                    @warn "Could not process $filepath during cleanup: $e"
                end
            end
        end
    end
    
    println("✅ Cleaned up $removed_count orphaned actor files.")
    return removed_count
end

"""
    create_knowledge_base(all_actors::Dict, categories::Dict, relationships::Dict) -> Dict

Create combined knowledge base from individual actor analyses with relationship data.
"""
function create_knowledge_base(all_actors::Dict, categories::Dict, relationships::Dict)
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
        ),
        "relationships" => relationships
    )
end

# Convenience function
"""
    extract_fuse_knowledge()

Extract knowledge using default paths (assumes running from knowledge/ directory).
"""
# Test function for a few actors
function test_extraction()
    extract_all_actors_parallel("../src", ".", limit=3, batch_size=3)
end

# Parallel extraction with default batch size
function extract_fuse_knowledge_parallel()
    extract_all_actors_parallel("../src", ".")
end

# Sequential extraction (for comparison)
function extract_fuse_knowledge()
    extract_all_actors("../src", ".")
end

# Utility functions for maintenance
function check_orphaned_fuse_actors()
    check_orphaned_actors("actors", "../src")
end

function clean_orphaned_fuse_actors()
    clean_orphaned_actors("actors", "../src")
end