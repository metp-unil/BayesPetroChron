using Turing
using Distributions
using MCMCChains
using Statistics
using BayesPetroChron


"""
    run_strict_ordering(obs_ages, obs_σ; T_min=50.0, T_max=300.0, Δ_max=200.0,
                         n_samples=3000, n_chains=4, target_accept=0.9)

Sample the strictly ordered age-sequence model:
grt_core (T1) ≥ grt_rim (T2) ≥ mnz (T3) ≥ zrn (T4) ≥ rt (T5)

"""
function run_strict_ordering(obs_ages, obs_σ; T_min=50.0, T_max=300.0, Δ_max=200.0,
                              n_samples=3000, n_chains=4, target_accept=0.9)

    @model function age_sequence(t_obs, σ_obs)
        T1  ~ Uniform(T_min, T_max)   # garnet core — oldest anchor
        Δ12 ~ Uniform(0.0, Δ_max)
        Δ23 ~ Uniform(0.0, Δ_max)
        Δ34 ~ Uniform(0.0, Δ_max)
        Δ45 ~ Uniform(0.0, Δ_max)

        T2 = T1 - Δ12   # garnet rim
        T3 = T2 - Δ23   # monazite
        T4 = T3 - Δ34   # zircon
        T5 = T4 - Δ45   # rutile — youngest
        T5 < T_min && (Turing.@addlogprob! -Inf; return)

        t_obs[1] ~ Normal(T1, σ_obs[1])   # grt_core
        t_obs[2] ~ Normal(T2, σ_obs[2])   # grt_rim
        t_obs[3] ~ Normal(T3, σ_obs[3])   # mnz
        t_obs[4] ~ Normal(T4, σ_obs[4])   # zrn
        t_obs[5] ~ Normal(T5, σ_obs[5])   # rt
    end

    model = age_sequence(obs_ages, obs_σ)
    chain = sample(model, NUTS(target_accept), MCMCThreads(), n_samples, n_chains)

    # Sample prior (oversample, 500k), then filter to keep only draws with T5 ≥ T_min
    chain_prior = sample(model, Prior(), 500_000)

    T1_prior  = vec(chain_prior[:T1])
    Δ12_prior = vec(chain_prior[:Δ12])
    Δ23_prior = vec(chain_prior[:Δ23])
    Δ34_prior = vec(chain_prior[:Δ34])
    Δ45_prior = vec(chain_prior[:Δ45])

    T2_prior = T1_prior .- Δ12_prior
    T3_prior = T2_prior .- Δ23_prior
    T4_prior = T3_prior .- Δ34_prior
    T5_prior = T4_prior .- Δ45_prior

    valid = T5_prior .>= T_min
    prior_vec = [T1_prior[valid], T2_prior[valid], T3_prior[valid], T4_prior[valid], T5_prior[valid]]

    # Recalculate posterior ages
    T1_samp  = vec(chain[:T1])
    Δ12_samp = vec(chain[:Δ12])
    Δ23_samp = vec(chain[:Δ23])
    Δ34_samp = vec(chain[:Δ34])
    Δ45_samp = vec(chain[:Δ45])

    T2_samp = T1_samp .- Δ12_samp
    T3_samp = T2_samp .- Δ23_samp
    T4_samp = T3_samp .- Δ34_samp
    T5_samp = T4_samp .- Δ45_samp

    post_samps = (grt_core = T1_samp, grt_rim = T2_samp, mnz = T3_samp, zrn = T4_samp, rt = T5_samp)

    dur_samps = (
        core_to_rim = T1_samp .- T2_samp,
        rim_to_mnz  = T2_samp .- T3_samp,
        mnz_to_zrn  = T3_samp .- T4_samp,
        zrn_to_rt   = T4_samp .- T5_samp,
        total       = T1_samp .- T5_samp,
    )

    post_vec = [T1_samp, T2_samp, T3_samp, T4_samp, T5_samp]

    # Chain-shaped (iteration × chain) counterparts of post_vec, needed for correct
    # ESS/R̂ on the derived ages (rim, mnz, zrn)
    post_mat = [
        Array(chain[:T1]),
        Array(chain[:T1] .- chain[:Δ12]),
        Array(chain[:T1] .- chain[:Δ12] .- chain[:Δ23]),
        Array(chain[:T1] .- chain[:Δ12] .- chain[:Δ23] .- chain[:Δ34]),
        Array(chain[:T1] .- chain[:Δ12] .- chain[:Δ23] .- chain[:Δ34] .- chain[:Δ45]),
    ]

    param_names  = [:T1, :Δ12, :Δ23, :Δ34, :Δ45]
    param_labels = ["T₁ (core)", "Δ₁₂", "Δ₂₃", "Δ₃₄", "Δ₄₅"]

    return (; chain, chain_prior, post_samps, dur_samps, post_vec, post_mat, prior_vec,
              param_names, param_labels, n_chains)
end

if abspath(PROGRAM_FILE) == @__FILE__

    include("geochronolgy_data.jl")

    results = cached_run("01_strict_ordering") do
        run_strict_ordering(obs_ages, obs_σ)
    end

    println("=== Convergence diagnostics ===")
    display(summarystats(results.chain))
    println()

    n_divergent = sum(results.chain[:numerical_error])
    n_total     = length(results.chain[:numerical_error])
    println("Divergent transitions: $n_divergent / $n_total " *
            "($(round(100 * n_divergent / n_total, digits=2))%)")
    println()

    println("=== Posterior age summaries (Ma) ===")
    for (name, s) in pairs(results.post_samps)
        lo, hi = quantile(s, [0.025, 0.975])
        println("  $name: median=$(round(median(s), digits=2))  " *
                "95% CI=[$(round(lo,digits=2)), $(round(hi,digits=2))]")
    end
    println()

    println("=== Posterior duration summaries (Ma) ===")
    for (name, s) in pairs(results.dur_samps)
        lo, hi = quantile(s, [0.025, 0.975])
        println("  $name: median=$(round(median(s), digits=2))  " *
                "95% CI=[$(round(lo,digits=2)), $(round(hi,digits=2))]")
    end

    fig = diagnostics_fig(results.chain; min_labels, obs_ages, obs_σ,
                           results.post_vec, results.post_mat, results.param_names,
                           results.param_labels, results.n_chains)
    display(fig)
end
