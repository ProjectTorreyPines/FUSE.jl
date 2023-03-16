import LibGit2
import Pkg
import TimeZones
import TimeZones: @tz_str
import Dates

"""
    warmup()

Function used to precompile the majority of FUSE
"""
function warmup()
    dd = IMAS.dd()
    return warmup(dd)
end

function warmup(dd::IMAS.dd)
    TimerOutputs.reset_timer!(to, "warmup")
    return TimerOutputs.@timeit to "warmup" begin
        TimerOutputs.@timeit to "init" begin
            ini, act = case_parameters(:FPP; version=:v1_demount, init_from=:scalars)
            init(dd, ini, act)
        end
        ActorWholeFacility(dd, act)
        return IMAS.freeze(dd)
    end
end

function get_remote_url(repo_path::AbstractString, remote_name::AbstractString="origin")
    repo = LibGit2.GitRepo(repo_path)

    # Get the remote URL
    remote = LibGit2.get(LibGit2.GitRemote, repo, remote_name)
    remote_url = LibGit2.url(remote)

    return remote_url
end


function recent_activity(repo_path::AbstractString; verbose=true)
    repo_name=basename(repo_path)
    repo = LibGit2.GitRepo(repo_path)

    # Check if the repository has any commits
    if LibGit2.isempty(repo)
        if verbose
            @info("$repo_name: The repository has no commits.")
        end
        return false
    end

    # Check if the repository is in a detached HEAD state
    try
        LibGit2.head(repo)
    catch
        if verbose
            @info("$repo_name: The repository is in a detached HEAD state.")
        end
        return false
    end

    # Get the current branch reference and its commit time
    branch_ref = LibGit2.head(repo)
    head_commit = LibGit2.peel(LibGit2.GitCommit, branch_ref)
    commit_time = LibGit2.committer(head_commit).time

    # Convert the commit time (in seconds since epoch) to a DateTime object
    commit_time_dt = Dates.unix2datetime(commit_time) # UTC

    # Get the current time and calculate the difference
    now = TimeZones.DateTime(Dates.now(), tz"UTC")
    time_diff = now - commit_time_dt

    if LibGit2.isdirty(repo)
        if verbose
            @info("$repo_name: The workspace is dirty")
        end
        return true
    end

    # Check if there was a recent commit
    if time_diff < Dates.Minute(10)
        if verbose
            @info("$repo_name: Latest commit is less 10 mins old ($(time_diff))")
        end
        return true
    else
        return false
    end
end

function smart_precompile()
    if VERSION <= v"1.8"
        return
    end

    dev_dir = "/Users/meneghini/Coding/julia"#dirname(dirname(@__DIR__))

    modified_packages = String[]
    for dir in map(p -> p.name, values(Pkg.dependencies()))
        repo_path = joinpath(dev_dir,dir)
        if isdir(joinpath(repo_path,".git"))
            if recent_activity(repo_path)
                if contains(get_remote_url(repo_path),"ProjectTorreyPines")
                    push!(modified_packages, dir)
                end
            end
        end
    end

    if isempty(modified_packages)
        FUSE.warmup()
    else
        @info("FUSE is in active development. skipping SnoopPrecompile!")
    end
end

SnoopPrecompile.@precompile_setup begin
    SnoopPrecompile.@precompile_all_calls begin
        smart_precompile()
    end
end