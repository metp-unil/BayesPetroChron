using BayesPetroChron

include(joinpath(@__DIR__, "..", "age_models", "02_soft_ordering.jl"))
include(joinpath(@__DIR__, "..", "age_models", "geochronolgy_data.jl"))

soft = cached_run("02_soft_ordering") do
    run_soft_ordering(obs_ages, obs_σ; σ_rim=σ_ages.grt_rim)
end

for (name, data) in pairs(soft.post_samps)
    println("***$(name)")
    println("5th percentile: $(quantile(data, 0.05))")
    println("95th percentile: $(quantile(data, 0.95))")
    println("")
end

Δgrt = soft.dur_samps.core_to_rim
println("***Δgrt")
println("5th percentile: $(quantile(Δgrt, 0.05))")
println("95th percentile: $(quantile(Δgrt, 0.95))")
