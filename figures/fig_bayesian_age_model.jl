using CairoMakie
using BayesPetroChron
using BayesPetroChron: normal_curve, kde_curve, add_ci_bars!

include(joinpath(@__DIR__, "..", "age_models", "02_soft_ordering.jl"))
include(joinpath(@__DIR__, "..", "age_models", "geochronolgy_data.jl"))

soft = cached_run("02_soft_ordering") do
    run_soft_ordering(obs_ages, obs_σ; σ_rim=σ_ages.grt_rim)
end

Δgrt = soft.dur_samps.core_to_rim

@model function geotherm(P_grtcore, T_grtcore, P_grtrim, T_grtrim)
    P_grt_core ~ Uniform(P_grtcore[1], P_grtcore[2])
    T_grt_core ~ Uniform(T_grtcore[1], T_grtcore[2])
    P_grt_rim  ~ Uniform(P_grtrim[1], P_grtrim[2])
    T_grt_rim  ~ Uniform(T_grtrim[1], T_grtrim[2])
end

P_grtcore_min, P_grtcore_max = 1.2  , 2.0         # GPa
T_grtcore_min, T_grtcore_max = 520. , 555.        # °C
P_grtrim_min, P_grtrim_max   = 2.0  , 2.5         # GPa
T_grtrim_min, T_grtrim_max   = 566. , 577.        # °C

model_geotherm = geotherm((P_grtcore_min, P_grtcore_max), (T_grtcore_min, T_grtcore_max),
                          (P_grtrim_min, P_grtrim_max), (T_grtrim_min, T_grtrim_max))
n_samples = 3000
n_chains = 4
chain_geotherm = sample(model_geotherm, NUTS(), MCMCThreads(), n_samples, n_chains)
summarystats(chain_geotherm)

P_core_samp = vec(chain_geotherm[:P_grt_core])
T_core_samp = vec(chain_geotherm[:T_grt_core])
P_rim_samp  = vec(chain_geotherm[:P_grt_rim])
T_rim_samp  = vec(chain_geotherm[:T_grt_rim])

dTdP_samp = (T_core_samp .- T_rim_samp) ./ (P_core_samp .- P_rim_samp)

const lithostatic_gradient = 0.031      # GPa/km

ΔP = P_core_samp .- P_rim_samp
burial_rate = ΔP ./ lithostatic_gradient ./ Δgrt .* 10.     # cmyr-1


slab_dip_μ = 30.0
slab_dip_σ = 5.0

slab_dip = randn(12000) * slab_dip_σ .+ slab_dip_μ
slab_dip = slab_dip .* π ./ 180.0

horizontal_plate_velocity = burial_rate .* tan.(slab_dip)    # cmyr-1

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

fig = Figure(size=(800, 600))

# Measurements
ax_meas = Axis(fig[1, 1], ygridvisible=false, ylabel = "Density [a.u.]")
hideydecorations!(ax_meas, label=false)
plot_phases!(ax_meas, [normal_curve(μ, σ) for (μ, σ) in zip(obs_ages, obs_σ)]; with_legend=true)
# add_ci_bars!(ax_meas, [rand(Normal(μ, σ), 10_000) for (μ, σ) in zip(obs_ages, obs_σ)], min_colors; y_pos=0.8)
for (μ, σ, color) in zip(obs_ages, obs_σ, min_colors)
    text!(ax_meas, μ, pdf(Normal(μ, σ), μ); text="$(round(μ, digits=1))\nI", align=(:center, :bottom), fontsize=11, color=color)
end
text!(ax_meas, 0.02, 0.98; text="Geochronological\nages", space=:relative, align=(:left, :top), fontsize=16, font=:bold)
ylims!(ax_meas, nothing, 0.8)

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
ylims!(ax_model, nothing, 0.8)

linkxaxes!(ax_meas, ax_model)


ax_Δgrt = Axis(fig[1, 2], ygridvisible=false, xlabel = "Δ Garnet age [Ma]", ylabel = "Density [a.u.]")
hideydecorations!(ax_Δgrt, label=false)
# hist!(ax_Δgrt, Δgrt; bins=100, label="Δ Garnet age")
xs, ys = kde_curve(Δgrt, boundary=(minimum(Δgrt), 2 * maximum(Δgrt)))
lines!(ax_Δgrt, xs, ys; color=:black, linewidth=2)
band!(ax_Δgrt, xs, zeros(length(xs)), ys; color=(:black, 0.2))
# add_ci_bars!(ax_Δgrt, [Δgrt], [:black]; y_pos=0.1)
let i = argmax(ys)
    text!(ax_Δgrt, xs[i], ys[i]; text="$(round(xs[i], digits=1))\nI", align=(:center, :bottom), fontsize=11, color=:black)
end
xlims!(ax_Δgrt, minimum(Δgrt) - 2, maximum(Δgrt) + 2)
ylims!(ax_Δgrt, nothing, 0.09)

ax_tectrate = Axis(fig[2, 2], ygridvisible=false, xlabel = "Tectonic rate [cmyr⁻¹]", ylabel = "Density [a.u.]")
hideydecorations!(ax_tectrate, label=false)
# hist!(ax_tectrate, burial_rate; bins=range(-100, 1, length=100), label="Tectonic rate", alpha=0.5)
# hist!(ax_tectrate, horizontal_plate_velocity; bins=range(-100, 1, length=100), label="Horizontal plate velocity", alpha=0.5)
xs_br, ys_br = kde_curve(burial_rate, boundary=(-200, maximum(burial_rate)))
xs_hpv, ys_hpv = kde_curve(horizontal_plate_velocity, boundary=(-200, maximum(horizontal_plate_velocity)))
lines!(ax_tectrate, xs_br, ys_br; color=:blue, linewidth=2, label="Burial rate")
band!(ax_tectrate, xs_br, zeros(length(xs_br)), ys_br; color=(:blue, 0.2))
let i = argmax(ys_br)
    text!(ax_tectrate, xs_br[i], ys_br[i]; text="$(round(xs_br[i], digits=2))\nI", align=(:center, :bottom), fontsize=11, color=:blue)
end
lines!(ax_tectrate, xs_hpv, ys_hpv; color=:red, linewidth=2, label="Horizontal plate velocity")
band!(ax_tectrate, xs_hpv, zeros(length(xs_hpv)), ys_hpv; color=(:red, 0.2))
let i = argmax(ys_hpv)
    text!(ax_tectrate, xs_hpv[i], ys_hpv[i]; text="$(round(xs_hpv[i], digits=2))\nI", align=(:center, :bottom), fontsize=11, color=:red)
end
# add_ci_bars!(ax_tectrate, [tectonic_rate, horizontal_plate_velocity], [:blue, :red]; y_pos=0.08, dodge_gap=0.005)

xlims!(ax_tectrate, -100, 4)
ylims!(ax_tectrate, nothing, 0.045)
ax_tectrate.xticks = LinearTicks(6)

# add legend
axislegend(ax_tectrate; position=:lt, framevisible=false, labelsize=10, orientation=:vertical)

fig
save(joinpath(@__DIR__, "..", "figures", "bayesian_age_model.pdf"), fig)