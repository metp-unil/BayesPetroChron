using CairoMakie
using BayesPetroChron
using BayesPetroChron: normal_curve, kde_curve

include(joinpath(@__DIR__, "..", "age_models", "02_soft_ordering.jl"))
include(joinpath(@__DIR__, "..", "age_models", "geochronolgy_data.jl"))

soft = cached_run("02_soft_ordering") do
    run_soft_ordering(obs_ages, obs_σ; σ_rim=σ_ages.grt_rim)
end

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

# HELPER: Ridgeline plot — each distribution normalized to its own peak and
# stacked vertically. `spacing` controls the vertical gap between baselines;
# values < 1.0 produce overlap. Phase names are shown as ytick labels.
function plot_ridgeline!(ax, curves, labels, colors; spacing=1.2, linewidth=2)
    baselines = [(i - 1) * spacing for i in eachindex(curves)]
    for ((xs, ys), color, baseline) in zip(curves, colors, baselines)
        ys_norm = ys ./ maximum(ys)
        band!(ax, xs, fill(baseline, length(xs)), ys_norm .+ baseline; color=(color, 0.3))
        lines!(ax, xs, ys_norm .+ baseline; color=color, linewidth=linewidth)
    end
    ax.yticks = (baselines, labels)
end

fig = Figure(size=(800, 600))

row_titles = ["Measured" "Soft ordering"]

# Row 1 — measurements
ax_meas = Axis(fig[1, 1], ygridvisible=false, ylabel = "Density [a.u.]")
hideydecorations!(ax_meas, label=false)
plot_phases!(ax_meas, [normal_curve(μ, σ) for (μ, σ) in zip(obs_ages, obs_σ)]; with_legend=true)
text!(ax_meas, 0.02, 0.98; text=row_titles[1], space=:relative, align=(:left, :top), fontsize=16, font=:bold)

# Row 1, col 2 — legend
Legend(fig[1, 2], ax_meas; framevisible=false, labelsize=16, labelfont=:bold, linewidth=6, tellwidth=false)

# Row 2 — one model: posterior (col 1), prior as ridgeline (col 2)
ax_post = Axis(fig[2, 1],
    ygridvisible=false,
    title  = "Posterior",
    xlabel = "Age [Ma]",
    ylabel = "Density [a.u.]")
hideydecorations!(ax_post, label=false)
plot_phases!(ax_post, [kde_curve(s) for s in soft.post_vec])
text!(ax_post, 0.02, 0.98; text=row_titles[2], space=:relative, align=(:left, :top), fontsize=16, font=:bold)

ax_prior = Axis(fig[2, 2],
    ygridvisible=false,
    title  = "Prior",
    xlabel = "Age [Ma]",
    ylabel = "Density [a.u.]")
hideydecorations!(ax_prior, label=false)
prior_curves = [kde_curve(s) for s in soft.prior_vec]
plot_ridgeline!(ax_prior, prior_curves, min_labels, min_colors)


linkxaxes!(ax_meas, ax_post)

fig
save(joinpath(@__DIR__, "..", "figures", "bayesian_age_model_with_priors.pdf"), fig)