const RESULTS_DIR = joinpath(@__DIR__, "..", "results")

"""
    cached_run(f, name; force=false)

Call zero-argument `f()` and cache its result to `results/<name>.jld2`. If a
cache file already exists and `force=false`, load and return the cached
result instead of recomputing (skipping, e.g., a full NUTS resampling run).
Pass `force=true` to recompute and overwrite the cache.

Designed for `do`-block use:
```julia
results = cached_run("01_strict_ordering") do
    run_strict_ordering(obs_ages, obs_σ)
end
```

`f` is expected to return a `NamedTuple` (e.g. an age-model's `run_*`
function), so downstream scripts can consume `results.chain`, `results.post_vec`,
etc. regardless of whether the value was freshly computed or loaded from disk.
"""
function cached_run(f, name; force=false)
    path = joinpath(RESULTS_DIR, "$name.jld2")
    if !force && isfile(path)
        @info "Loading cached results for \"$name\" from $path"
        return load(path, "results")
    end
    results = f()
    mkpath(RESULTS_DIR)
    @info "Saving results for \"$name\" to $path"
    save(path, "results", results)
    return results
end
