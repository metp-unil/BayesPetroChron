using Turing
using Distributions
using MCMCChains
using Statistics
using BayesPetroChron


"""
    run_soft_ordering(obs_ages, obs_σ; σ_rim, T_min=50.0, T_max=300.0, Δ_max=200.0,
                       n_samples=3000, n_chains=4, target_accept=0.95)

Sample the soft-ordered age-sequence model:
grt_core (T1) ≥ grt_rim (T2)
grt_core (T1) ≥ mnz (T3)
grt_core (T1) ≥ zrn (T4)
rt (T5) — unconstrained


Grt rim, mnz, and zrn may have grown during the same metamorphic event: their relative order is
not hard-enforced, but a reversal (mnz or zrn older than grt rim) is penalised in units of σ_rim.
The penalty is quadratic, so a reversal within ~1σ costs almost nothing, while one of several σ
is strongly disfavoured.
"""
function run_soft_ordering(obs_ages, obs_σ; σ_rim, T_min=50.0, T_max=300.0, Δ_max=200.0,
                            n_samples=3000, n_chains=4, target_accept=0.95)

    @model function age_sequence_soft(t_obs, σ_obs)
        T1  ~ Uniform(T_min, T_max)   # garnet core
        T5  ~ Uniform(T_min, T_max)   # rutile

        Δ12 ~ Uniform(0.0, Δ_max)
        T2 = T1 - Δ12   # garnet rim

        Δ13 ~ Uniform(0.0, Δ_max)
        Δ14 ~ Uniform(0.0, Δ_max)
        T3 = T1 - Δ13   # monazite
        T4 = T1 - Δ14   # zircon

        # check Ti < T_min, otherwise reject sample
        any(T -> T < T_min, (T2, T3, T4)) && (Turing.@addlogprob! -Inf; return)

        # Penalise mnz > grt_rim / zrn > grt_rim, quadratically in units of σ_rim — no hard
        # rejection, just a smooth preference.
        reversal_Mnz = max(0.0, T3 - T2)
        reversal_Zrn = max(0.0, T4 - T2)
        Turing.@addlogprob! -(reversal_Mnz / σ_rim)^2
        Turing.@addlogprob! -(reversal_Zrn / σ_rim)^2

        t_obs[1] ~ Normal(T1, σ_obs[1])   # grt_core
        t_obs[2] ~ Normal(T2, σ_obs[2])   # grt_rim
        t_obs[3] ~ Normal(T3, σ_obs[3])   # mnz
        t_obs[4] ~ Normal(T4, σ_obs[4])   # zrn
        t_obs[5] ~ Normal(T5, σ_obs[5])   # rt
    end

    model = age_sequence_soft(obs_ages, obs_σ)
    chain = sample(model, NUTS(target_accept), MCMCThreads(), n_samples, n_chains)

    # Sample prior (oversample, 500k), then filter to keep only draws with T5 ≥ T_min
    # add soft reversal penalty via accept/reject step with probability exp(-(reversal/σ_rim)^2)
    chain_prior = sample(model, Prior(), 500_000)

    T1_prior  = vec(chain_prior[:T1])
    T5_prior  = vec(chain_prior[:T5])
    Δ12_prior = vec(chain_prior[:Δ12])
    Δ13_prior = vec(chain_prior[:Δ13])
    Δ14_prior = vec(chain_prior[:Δ14])

    T2_prior = T1_prior .- Δ12_prior
    T3_prior = T1_prior .- Δ13_prior
    T4_prior = T1_prior .- Δ14_prior

    hard_ok = (T2_prior .>= T_min) .& (T3_prior .>= T_min) .& (T4_prior .>= T_min)
    reversal_mnz_prior = max.(0.0, T3_prior .- T2_prior)
    reversal_zrn_prior = max.(0.0, T4_prior .- T2_prior)
    accept_prob = exp.(-(reversal_mnz_prior ./ σ_rim) .^ 2) .* exp.(-(reversal_zrn_prior ./ σ_rim) .^ 2)
    valid = hard_ok .& (rand(length(T1_prior)) .< accept_prob)

    prior_vec = [T1_prior[valid], T2_prior[valid], T3_prior[valid], T4_prior[valid], T5_prior[valid]]

    # Recalculate posterior ages
    T1_samp  = vec(chain[:T1])
    Δ12_samp = vec(chain[:Δ12])
    Δ13_samp = vec(chain[:Δ13])
    Δ14_samp = vec(chain[:Δ14])
    T5_samp  = vec(chain[:T5])

    T2_samp = T1_samp .- Δ12_samp
    T3_samp = T1_samp .- Δ13_samp
    T4_samp = T1_samp .- Δ14_samp

    post_samps = (grt_core = T1_samp, grt_rim = T2_samp, mnz = T3_samp, zrn = T4_samp, rt = T5_samp)

    dur_samps = (
        core_to_rim = T1_samp .- T2_samp,
        rim_to_mnz  = T2_samp .- T3_samp,
        rim_to_zrn  = T2_samp .- T4_samp,
        mnz_to_rt   = T3_samp .- T5_samp,
        zrn_to_rt   = T4_samp .- T5_samp,
        total       = T1_samp .- T5_samp,
    )

    p_mnz_older_rim  = mean(dur_samps.rim_to_mnz .< 0)
    p_mnz_older_core = mean(T1_samp .- T3_samp .< 0)
    p_zrn_older_rim  = mean(dur_samps.rim_to_zrn .< 0)
    p_zrn_older_core = mean(T1_samp .- T4_samp .< 0)

    post_vec = [T1_samp, T2_samp, T3_samp, T4_samp, T5_samp]

    # Chain-shaped (iteration × chain) counterparts of post_vec, needed for correct
    # ESS/R̂ on the derived ages (rim, mnz, zrn)
    post_mat = [
        Array(chain[:T1]),
        Array(chain[:T1] .- chain[:Δ12]),
        Array(chain[:T1] .- chain[:Δ13]),
        Array(chain[:T1] .- chain[:Δ14]),
        Array(chain[:T5]),
    ]

    param_names  = [:T1, :Δ12, :Δ13, :Δ14, :T5]
    param_labels = ["T₁ (core)", "Δ₁₂ (core→rim)", "Δ₁₃ (core→mnz)", "Δ₁₄ (core→zrn)", "T₅ (rt)"]

    return (; chain, chain_prior, post_samps, dur_samps, post_vec, post_mat, prior_vec,
              param_names, param_labels, n_chains,
              p_mnz_older_rim, p_mnz_older_core, p_zrn_older_rim, p_zrn_older_core)
end

if abspath(PROGRAM_FILE) == @__FILE__
    include("geochronolgy_data.jl")

    results = cached_run("02_soft_ordering") do
        run_soft_ordering(obs_ages, obs_σ; σ_rim=σ_ages.grt_rim)
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
    println()

    println("P(mnz older than grt_rim)  ≈ $(round(results.p_mnz_older_rim, digits=3))")
    println("P(mnz older than grt_core) ≈ $(round(results.p_mnz_older_core, digits=3))")
    println("P(zrn older than grt_rim)  ≈ $(round(results.p_zrn_older_rim, digits=3))")
    println("P(zrn older than grt_core) ≈ $(round(results.p_zrn_older_core, digits=3))")

    fig = diagnostics_fig(results.chain; min_labels, obs_ages, obs_σ,
                           results.post_vec, results.post_mat, results.param_names,
                           results.param_labels, results.n_chains)
    display(fig)
end
