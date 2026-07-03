using CairoMakie
using BayesPetroChron
using BayesPetroChron: normal_curve, kde_curve

include(joinpath(@__DIR__, "..", "age_models", "01_strict_ordering.jl"))
include(joinpath(@__DIR__, "..", "age_models", "02_soft_ordering.jl"))
include(joinpath(@__DIR__, "..", "age_models", "geochronolgy_data.jl"))

strict = cached_run("01_strict_ordering") do
    run_strict_ordering(obs_ages, obs_σ)
end

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

fig = Figure(size=(800, 800))

row_titles = ["Measured", "Strict ordering", "Soft ordering"]

row_axes = [Axis[] for _ in 1:3]

# Row 1 — measurements
ax_meas = Axis(fig[1, 1], ygridvisible=false, ylabel = "Density [a.u.]")
hideydecorations!(ax_meas, label=false)
plot_phases!(ax_meas, [normal_curve(μ, σ) for (μ, σ) in zip(obs_ages, obs_σ)]; with_legend=true)
text!(ax_meas, 0.02, 0.98; text=row_titles[1], space=:relative, align=(:left, :top), fontsize=16, font=:bold)
push!(row_axes[1], ax_meas)

# Row 1, col 2 — legend
Legend(fig[1, 2], ax_meas; framevisible=false, labelsize=16, labelfont=:bold, linewidth=6, tellwidth=false)

# Rows 2–3 — one model per row: posterior (col 1), prior as ridgeline (col 2)
for (row, results) in zip(2:3, (strict, soft))
    is_last_row = row == 3

    ax_post = Axis(fig[row, 1],
        ygridvisible=false,
        title  = row == 2 ? "Posterior" : "",
        xlabel = is_last_row ? "Age [Ma]" : "",
        ylabel = "Density [a.u.]")
    hideydecorations!(ax_post, label=false)
    plot_phases!(ax_post, [kde_curve(s) for s in results.post_vec])
    text!(ax_post, 0.02, 0.98; text=row_titles[row], space=:relative, align=(:left, :top), fontsize=16, font=:bold)
    push!(row_axes[row], ax_post)

    ax_prior = Axis(fig[row, 2],
        ygridvisible=false,
        title  = row == 2 ? "Prior" : "",
        xlabel = is_last_row ? "Age [Ma]" : "",
        ylabel = "Density [a.u.]")
    hideydecorations!(ax_prior, label=false)
    prior_curves = [kde_curve(s) for s in results.prior_vec]
    plot_ridgeline!(ax_prior, prior_curves, min_labels, min_colors)
    push!(row_axes[row], ax_prior)
end

linkxaxes!(row_axes[1][1], row_axes[2][1], row_axes[3][1])
linkxaxes!(row_axes[2][2], row_axes[3][2])

fig
save(joinpath(@__DIR__, "..", "figures", "bayesian_age_model_ridgeline.pdf"), fig)