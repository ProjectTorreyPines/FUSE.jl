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
    generate_github_url(actor_file::String, src_dir::String) -> String

Generate a GitHub URL for the given actor source file.
Assumes the FUSE.jl repository structure and master branch.
"""
function generate_github_url(actor_file::String, src_dir::String)
    # Calculate relative path from src directory
    rel_path = relpath(actor_file, src_dir)
    
    # Base GitHub URL for FUSE.jl repository
    base_url = "https://github.com/ProjectTorreyPines/FUSE.jl/blob/master/src"
    
    # Combine to create full URL
    return joinpath(base_url, rel_path)
end

"""
    extract_all_actors(src_dir::String, knowledge_dir::String="knowledge")

Extract knowledge for all FUSE actors using Claude analysis.
Creates individual JSON files and combined knowledge base.
"""
function extract_all_actors_parallel(src_dir::String, knowledge_dir::String="knowledge"; batch_size=8, limit=nothing, skip_existing=true)
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
    
    # Find all actor files (optionally skip those with existing JSON files)
    actor_files = find_actor_files(actors_src_dir, actors_kb_dir, skip_existing)
    if limit !== nothing
        actor_files = actor_files[1:min(limit, length(actor_files))]
    end
    
    if skip_existing && length(actor_files) == 0
        println("✅ All actor JSON files are up to date! No extraction needed.")
        
        # Still generate combined knowledge base from existing files
        all_actors = Dict()
        categories = Dict()
        
        # Post-processing: Add GitHub URLs to any existing JSON files that don't have them
        add_source_urls_to_existing_jsons!(src_dir, knowledge_dir)
        
        missing_actors = include_existing_actors!(all_actors, categories, actors_kb_dir)
        relationships = analyze_actor_relationships(all_actors)
        knowledge_base = create_knowledge_base(all_actors, categories, relationships)
        
        kb_path = joinpath(knowledge_dir, "fuse_knowledge_base.json")
        open(kb_path, "w") do io
            JSON.print(io, knowledge_base, 2)
        end
        
        println("Updated combined knowledge base: $kb_path")
        println("Processed $(length(all_actors)) actors across $(length(categories)) categories")
        return knowledge_base
    end
    
    println("Found $(length(actor_files)) actor files$(skip_existing ? " needing extraction" : "")")
    
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
                    analyze_actor_with_claude_retry(actor_file, category, schema_path, src_dir)
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
    
    # Post-processing: Add GitHub URLs to any existing JSON files that don't have them
    url_updates = add_source_urls_to_existing_jsons!(src_dir, knowledge_dir)
    
    # Check for and include any existing actors from individual files not processed in this run
    missing_actors = include_existing_actors!(all_actors, categories, actors_kb_dir)
    
    # Check for orphaned individual files (JSONs without corresponding source files)
    orphaned_actors = check_orphaned_actors(actors_kb_dir, src_dir)
    
    if !isempty(missing_actors) || !isempty(orphaned_actors) || url_updates > 0
        if !isempty(missing_actors)
            println("📂 Included $(length(missing_actors)) existing actors in knowledge base: $(join(missing_actors, ", "))")
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
function extract_all_actors(src_dir::String, knowledge_dir::String="knowledge"; limit=nothing, skip_existing=true)
    return extract_all_actors_parallel(src_dir, knowledge_dir; batch_size=1, limit=limit, skip_existing=skip_existing)
end

"""
    find_actor_files(actors_dir::String, knowledge_dir::String="", skip_existing::Bool=false) -> Vector{Tuple{String, String}}

Find all *_actor.jl files and determine their categories from directory structure.
If skip_existing=true and knowledge_dir is provided, filters out files that already have JSON files.
Returns vector of (filepath, category) tuples.
"""
function find_actor_files(actors_dir::String, knowledge_dir::String="", skip_existing::Bool=false)
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
                
                # Check if we should skip files with existing JSON files
                if skip_existing && !isempty(knowledge_dir)
                    # Calculate expected JSON file path
                    rel_path = relpath(filepath, actors_dir)
                    json_filename = replace(basename(rel_path), ".jl" => ".json")
                    json_rel_path = joinpath(dirname(rel_path), json_filename)
                    json_filepath = joinpath(knowledge_dir, json_rel_path)
                    
                    # Skip if JSON file already exists
                    if isfile(json_filepath)
                        continue
                    end
                end
                
                push!(actor_files, (filepath, category))
            end
        end
    end
    
    return actor_files
end

"""
    analyze_actor_with_claude_retry(actor_file::String, category::String, schema_path::String, src_dir::String) -> Dict

Use Claude to analyze a single actor file with retry logic.
"""
function analyze_actor_with_claude_retry(actor_file::String, category::String, schema_path::String, src_dir::String; max_retries=1)
    for attempt in 1:(max_retries + 1)
        result = analyze_actor_with_claude(actor_file, category, schema_path, src_dir)
        
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
function analyze_actor_with_claude(actor_file::String, category::String, schema_path::String, src_dir::String="../src")
    # Read actor source code
    if !isfile(actor_file)
        return Dict("error" => "File not found: $actor_file")
    end
    
    source_code = read(actor_file, String)
    
    # Read schema for Claude
    schema_content = read(schema_path, String)
    
    # Create analysis prompt
    prompt = """
You are analyzing a FUSE actor

Please analyze this Julia source code and return a JSON object that exactly matches the provided schema.

NOTE: You can use the docstring for a high-level description of what the actor does,
but you should not trust that. The only source of truth is the Julia executable code.
Please always cross-check that the description of the actors, parameters and such
that you find ind the docstrings is indeed up to date. If any discrepancy comes up,
always trust the executable julia code.

**CRITICAL INSTRUCTIONS:**
- Return ONLY valid JSON, no markdown, no explanations, no code blocks
- The JSON must validate against the schema
- Fill out ALL required fields
- Use the category "$category" for the category field
- Extract the actual actor name from the struct definition in the code
- For sub_actors: only include if this is a compound actor (ie. inherits from CompoundAbstractActor)
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
        
        # Add GitHub source URL
        actor_data["source_url"] = generate_github_url(actor_file, src_dir)
        
        # Validate required fields exist (including the new source_url)
        required_fields = ["name", "category", "hierarchy", "description", "physics_domain", 
                          "data_inputs", "data_outputs", "sub_actors", "key_parameters", "usage_notes", "source_url"]
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
    include_existing_actors!(all_actors::Dict, categories::Dict, actors_kb_dir::String) -> Vector{String}

Check for actors that exist as individual JSON files but weren't included in the current extraction run.
Add any such actors to the combined knowledge base and return their names.
This ensures the combined knowledge base includes all available actor data, not just newly extracted ones.
"""
function include_existing_actors!(all_actors::Dict, categories::Dict, actors_kb_dir::String)
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
                        @info "Including existing actor in knowledge base: $actor_name from $filepath"
                        
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
    extract_all_actors_parallel("../src", ".", limit=3, batch_size=3, skip_existing=true)
end

# Parallel extraction with default batch size (skips existing by default)
function extract_fuse_knowledge_parallel()
    extract_all_actors_parallel("../src", ".", skip_existing=true)
end

# Sequential extraction (skips existing by default)
function extract_fuse_knowledge()
    extract_all_actors("../src", ".", skip_existing=true)
end

# Force re-extraction of all actors (override skip_existing)
function extract_fuse_knowledge_force()
    extract_all_actors_parallel("../src", ".", skip_existing=false)
end

# Utility functions for maintenance
function check_orphaned_fuse_actors()
    check_orphaned_actors("actors", "../src")
end

function clean_orphaned_fuse_actors()
    clean_orphaned_actors("actors", "../src")
end

"""
    add_source_urls_to_existing_jsons!(src_dir::String, knowledge_dir::String) -> Int

Post-processing step to add GitHub source URLs to existing JSON files that don't have them.
Returns the number of files updated.
"""
function add_source_urls_to_existing_jsons!(src_dir::String, knowledge_dir::String)
    actors_src_dir = joinpath(src_dir, "actors")
    actors_kb_dir = joinpath(knowledge_dir, "actors")
    
    if !isdir(actors_kb_dir)
        println("No knowledge directory found at: $actors_kb_dir")
        return 0
    end
    
    println("🔗 Adding GitHub source URLs to existing JSON files...")
    updated_count = 0
    
    # Walk through all existing JSON files
    for (root, dirs, files) in walkdir(actors_kb_dir)
        for file in files
            if endswith(file, ".json")
                json_filepath = joinpath(root, file)
                
                try
                    # Read existing JSON
                    actor_data = JSON.parsefile(json_filepath)
                    
                    # Check if source_url field is missing or empty
                    if !haskey(actor_data, "source_url") || isempty(actor_data["source_url"])
                        # Find corresponding source file
                        source_filename = replace(file, ".json" => ".jl")
                        rel_path = relpath(json_filepath, actors_kb_dir)
                        source_rel_path = joinpath(dirname(rel_path), source_filename)
                        source_filepath = joinpath(actors_src_dir, source_rel_path)
                        
                        if isfile(source_filepath)
                            # Generate GitHub URL
                            github_url = generate_github_url(source_filepath, src_dir)
                            
                            # Add source_url field
                            actor_data["source_url"] = github_url
                            
                            # Write updated JSON back to file
                            open(json_filepath, "w") do io
                                JSON.print(io, actor_data, 2)
                            end
                            
                            updated_count += 1
                            println("  ✓ Added URL to $(actor_data["name"]): $github_url")
                        else
                            @warn "Source file not found for JSON: $json_filepath (expected: $source_filepath)"
                        end
                    end
                    
                catch e
                    @warn "Failed to process JSON file $json_filepath: $e"
                end
            end
        end
    end
    
    println("📂 Updated $updated_count JSON files with GitHub URLs")
    return updated_count
end

# Convenience function for post-processing
function add_fuse_source_urls()
    add_source_urls_to_existing_jsons!("../src", ".")
end