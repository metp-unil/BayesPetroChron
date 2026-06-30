using Turing
using Distributions
using MCMCChains
using Statistics
using BayesPetroChron

# Measured ages (Ma)
μ_ages = (
    grt_core = 136.7,
    grt_rim  = 139.7,
    mnz      = 138.8,
    zrn      = 137.2,
    rt       = 121.5
)

σ_ages = (
    grt_core = 20.1,
    grt_rim  = 8.7,
    mnz      = 3.4,
    zrn      = 1.2,
    rt       = 2.5
)

obs_ages = [μ_ages.grt_core, μ_ages.grt_rim, μ_ages.mnz, μ_ages.zrn, μ_ages.rt]
obs_σ    = [σ_ages.grt_core, σ_ages.grt_rim, σ_ages.mnz, σ_ages.zrn, σ_ages.rt] .* 0.5    # reported as 2σ

labels = ["Garnet core", "Garnet rim", "Monazite", "Zircon", "Rutile"]
colors = [:firebrick, :orange, :seagreen, :steelblue, :purple]

# Age-sequence model
# grt_core (T1) ≥ grt_rim (T2) ≥ mnz (T3) ≥ zrn (T4) ≥ rt (T5)

const T_min = 50.0    # Ma — broad geological lower bound
const T_max = 300.0   # Ma — broad geological upper bound
const Δ_max = 200.0   # Ma — maximum plausible interval between events

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


n_samples = 3000
n_chains  = 4

chain = sample(model, NUTS(0.9), MCMCThreads(), n_samples, n_chains)

println("=== Convergence diagnostics ===")
display(summarystats(chain))
println()

n_divergent = sum(chain[:numerical_error])
n_total     = length(chain[:numerical_error])
println("Divergent transitions: $n_divergent / $n_total " *
        "($(round(100 * n_divergent / n_total, digits=2))%)")
println()

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

println("=== Posterior age summaries (Ma) ===")
for (name, s) in pairs(post_samps)
    lo, hi = quantile(s, [0.025, 0.975])
    println("  $name: median=$(round(median(s), digits=2))  " *
            "95% CI=[$(round(lo,digits=2)), $(round(hi,digits=2))]")
end
println()

dur_samps = (
    core_to_rim = T1_samp .- T2_samp,
    rim_to_mnz  = T2_samp .- T3_samp,
    mnz_to_zrn  = T3_samp .- T4_samp,
    zrn_to_rt   = T4_samp .- T5_samp,
    total       = T1_samp .- T5_samp,
)

println("=== Posterior duration summaries (Ma) ===")
for (name, s) in pairs(dur_samps)
    lo, hi = quantile(s, [0.025, 0.975])
    println("  $name: median=$(round(median(s), digits=2))  " *
            "95% CI=[$(round(lo,digits=2)), $(round(hi,digits=2))]")
end

# Plotting
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

fig = diagnostics_fig(chain; labels, obs_ages, obs_σ, post_vec, post_mat,
                       param_names, param_labels, n_chains)

fig_compare = meas_v_model_age_fig(; labels, colors, obs_ages, obs_σ, post_vec)

