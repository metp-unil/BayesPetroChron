using CairoMakie
using BayesPetroChron
using BayesPetroChron: normal_curve, kde_curve

include(joinpath(@__DIR__, "..", "age_models", "02_soft_ordering.jl"))
include(joinpath(@__DIR__, "..", "age_models", "geochronolgy_data.jl"))

soft_RG18 = cached_run("02_soft_ordering") do
    run_soft_ordering(obs_ages, obs_σ; σ_rim=σ_ages.grt_rim)
end

μ_ages_RG03 = (
    grt_core = 128.2,
    grt_rim  = 140.6,
    rt       = 125.6,
)

σ_ages_RG03 = (
    grt_core = 10.7,
    grt_rim  = 8.4,
    rt       = 2.7,
)

obs_ages_RG03 = [μ_ages_RG03.grt_core, μ_ages_RG03.grt_rim, μ_ages_RG03.rt]
obs_σ_RG03    = [σ_ages_RG03.grt_core, σ_ages_RG03.grt_rim, σ_ages_RG03.rt] .* 0.5    # reported as 2σ

@model function age_sequence_soft(t_obs, σ_obs)
    T_min=50.0
    T_max=250.0
    Δ_max=100.0

    T1  ~ Uniform(T_min, T_max)   # garnet core
    T3  ~ Uniform(T_min, T_max)   # rutile

    Δ12 ~ Uniform(0.0, Δ_max)
    T2 = T1 - Δ12   # garnet rim

    # check Ti < T_min, otherwise reject sample
    any(T -> T < T_min, (T2,)) && (Turing.@addlogprob! -Inf; return)

    t_obs[1] ~ Normal(T1, σ_obs[1])   # grt_core
    t_obs[2] ~ Normal(T2, σ_obs[2])   # grt_rim
    t_obs[3] ~ Normal(T3, σ_obs[3])   # rt
end

model_RG03 = age_sequence_soft(obs_ages_RG03, obs_σ_RG03)

n_samples=3000
n_chains=4

chain_RG03 = sample(model_RG03, NUTS(), MCMCThreads(), n_samples, n_chains)

summarystats(chain_RG03)

# Recalculate posterior ages
T1_samp_RG03  = vec(chain_RG03[:T1])
T2_samp_RG03  = T1_samp_RG03 .- vec(chain_RG03[:Δ12])
T3_samp_RG03  = vec(chain_RG03[:T3])






rg03_idx    = [1, 2, 5]   # core, rim, rutile indices into min_labels/min_colors
rg03_colors = min_colors[rg03_idx]
rg03_labels = min_labels[rg03_idx]
rg03_post_vec = [T1_samp_RG03, T2_samp_RG03, T3_samp_RG03]

fig = Figure(size=(1000, 600))

Label(fig[0, 1], "RG18"; fontsize=14, font=:bold, tellwidth=false)
Label(fig[0, 2], "RG03"; fontsize=14, font=:bold, tellwidth=false)

# ── Col 1: RG18 ───────────────────────────────────────────────────────────────
ax_meas = Axis(fig[1, 1], ygridvisible=false, ylabel = "Density [a.u.]")
hideydecorations!(ax_meas, label=false)
plot_phases!(ax_meas, [normal_curve(μ, σ) for (μ, σ) in zip(obs_ages, obs_σ)]; with_legend=true)
for (μ, σ, color) in zip(obs_ages, obs_σ, min_colors)
    text!(ax_meas, μ, pdf(Normal(μ, σ), μ); text="$(round(μ, digits=1))\nI", align=(:center, :bottom), fontsize=11, font=:bold, color=color)
end
text!(ax_meas, 0.02, 0.98; text="Geochronological\nages", space=:relative, align=(:left, :top), fontsize=16, font=:bold)
ylims!(ax_meas, nothing, 0.8)

ax_model = Axis(fig[2, 1], ygridvisible=false, xlabel = "Age [Ma]", ylabel = "Density [a.u.]")
hideydecorations!(ax_model, label=false)
plot_phases!(ax_model, [kde_curve(s) for s in soft_RG18.post_vec])
for (samps, color) in zip(soft_RG18.post_vec, min_colors)
    xs_m, ys_m = kde_curve(samps)
    i = argmax(ys_m)
    text!(ax_model, xs_m[i], ys_m[i]; text="$(round(xs_m[i], digits=1))\nI", align=(:center, :bottom), fontsize=11, font=:bold, color=color)
end
text!(ax_model, 0.02, 0.98; text="Inferred posterior\nages", space=:relative, align=(:left, :top), fontsize=16, font=:bold)
ylims!(ax_model, nothing, 0.8)

# ── Col 2: RG03 (grt core, grt rim, rutile only) ─────────────────────────────
ax_meas_RG03 = Axis(fig[1, 2], ygridvisible=false)
hideydecorations!(ax_meas_RG03)
for (μ, σ, color, label) in zip(obs_ages_RG03, obs_σ_RG03, rg03_colors, rg03_labels)
    xs_n, ys_n = normal_curve(μ, σ)
    lines!(ax_meas_RG03, xs_n, ys_n; color=color, linewidth=2)
    band!(ax_meas_RG03, xs_n, zeros(length(xs_n)), ys_n; color=(color, 0.2))
    text!(ax_meas_RG03, μ, pdf(Normal(μ, σ), μ); text="$(round(μ, digits=1))\nI", align=(:center, :bottom), fontsize=11, font=:bold, color=color)
end
text!(ax_meas_RG03, 0.02, 0.98; text="Geochronological\nages", space=:relative, align=(:left, :top), fontsize=16, font=:bold)
ylims!(ax_meas_RG03, nothing, 0.8)

ax_model_RG03 = Axis(fig[2, 2], ygridvisible=false, xlabel = "Age [Ma]")
hideydecorations!(ax_model_RG03)
for (samps, color) in zip(rg03_post_vec, rg03_colors)
    xs_m, ys_m = kde_curve(samps)
    lines!(ax_model_RG03, xs_m, ys_m; color=color, linewidth=2)
    band!(ax_model_RG03, xs_m, zeros(length(xs_m)), ys_m; color=(color, 0.2))
    i = argmax(ys_m)
    text!(ax_model_RG03, xs_m[i], ys_m[i]; text="$(round(xs_m[i], digits=1))\nI", align=(:center, :bottom), fontsize=11, font=:bold, color=color)
end
text!(ax_model_RG03, 0.02, 0.98; text="Inferred posterior\nages", space=:relative, align=(:left, :top), fontsize=16, font=:bold)
ylims!(ax_model_RG03, nothing, 0.8)

linkxaxes!(ax_meas, ax_model, ax_meas_RG03, ax_model_RG03)
linkyaxes!(ax_meas, ax_meas_RG03)
linkyaxes!(ax_model, ax_model_RG03)

fig

save(joinpath(@__DIR__, "..", "figures", "compare_samples.pdf"), fig)