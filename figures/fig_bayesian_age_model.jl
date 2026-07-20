using CairoMakie
using BayesPetroChron
using BayesPetroChron: normal_curve, kde_curve

include(joinpath(@__DIR__, "..", "age_models", "02_soft_ordering.jl"))
include(joinpath(@__DIR__, "..", "age_models", "geochronolgy_data.jl"))

soft = cached_run("02_soft_ordering") do
    run_soft_ordering(obs_ages, obs_σ; σ_rim=σ_ages.grt_rim)
end

print_results(soft, min_labels, to_file="age_model_results.txt")

# --------
# PLOTTING
# --------

# HELPER: Filled density curve per phase
function plot_phases!(ax, curves; linewidth=2, with_legend=false)
    for ((xs, ys), label, color) in zip(curves, min_labels, min_colors)
        lines!(ax, xs, ys; color=color, linewidth=linewidth, label = with_legend ? label : nothing)
        band!(ax, xs, zeros(length(xs)), ys; color=(color, 0.2))
    end
end

fig = Figure(size=(400, 600))

# Measurements
ax_meas = Axis(fig[1, 1], ygridvisible=false, ylabel = "Density [a.u.]")
hideydecorations!(ax_meas, label=false)
plot_phases!(ax_meas, [normal_curve(μ, σ) for (μ, σ) in zip(obs_ages, obs_σ)]; with_legend=true)
# add_ci_bars!(ax_meas, [rand(Normal(μ, σ), 10_000) for (μ, σ) in zip(obs_ages, obs_σ)], min_colors; y_pos=0.8)
for (μ, σ, color) in zip(obs_ages, obs_σ, min_colors)
    text!(ax_meas, μ, pdf(Normal(μ, σ), μ); text="$(round(μ, digits=1))\nI", align=(:center, :bottom), fontsize=11, color=color)
end
text!(ax_meas, 0.02, 0.98; text="Geochronological\nages", space=:relative, align=(:left, :top), fontsize=16, font=:bold)
ylims!(ax_meas, nothing, 0.45)

ax_model = Axis(fig[2, 1], ygridvisible=false, xlabel = "Age [Ma]", ylabel = "Density [a.u.]")
hideydecorations!(ax_model, label=false)
plot_phases!(ax_model, [kde_curve(s) for s in soft.post_vec])
# add_ci_bars!(ax_model, soft.post_vec, min_colors; y_pos=0.8)
for (samps, color) in zip(soft.post_vec, min_colors)
    xs_m, ys_m = kde_curve(samps)
    i = argmax(ys_m)
    text!(ax_model, xs_m[i], ys_m[i]; text="$(round(xs_m[i], digits=1))\nI", align=(:center, :bottom), fontsize=11, color=color)
end
text!(ax_model, 0.02, 0.98; text="Inferred posterior\nages", space=:relative, align=(:left, :top), fontsize=16, font=:bold)
ylims!(ax_model, nothing, 0.45)

linkxaxes!(ax_meas, ax_model)

fig
save(joinpath(@__DIR__, "..", "figures", "bayesian_age_model.pdf"), fig)