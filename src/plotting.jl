const FIGURE_DPI = 72

"""
    print_results(results, min_labels; to_file=nothing)

Print a formatted summary table of posterior age statistics for each mineral phase.
Columns: Median, Mode (KDE peak), MAD (median absolute deviation), 95% CI.

Pass a file path to `to_file` to additionally write the output to a text file.
"""
function print_results(results, min_labels; to_file=nothing)
    header = rpad("Mineral", 16) *
             lpad("Median", 8)   *
             lpad("Mode",   8)   *
             lpad("MAD",    8)   *
             lpad("95% CI", 22)
    rule   = "─" ^ length(header)

    rows = String[]
    for (label, samps) in zip(min_labels, results.post_vec)
        med      = median(samps)
        xs, ys   = kde_curve(samps)
        mode_val = xs[argmax(ys)]
        mad_val  = median(abs.(samps .- med))
        lo, hi   = quantile(samps, [0.025, 0.975])
        ci_str   = "[$(round(lo, digits=1)), $(round(hi, digits=1))]"
        push!(rows,
            rpad(label, 16)                          *
            lpad(string(round(med,      digits=1)), 8) *
            lpad(string(round(mode_val, digits=1)), 8) *
            lpad(string(round(mad_val,  digits=1)), 8) *
            lpad(ci_str, 22))
    end

    title  = "Posterior age summaries [Ma]"
    title_centered = lpad(title, div(length(header) + length(title), 2))
    output = join([title_centered, "", header, rule, rows..., ""], "\n")
    print(output)

    if !isnothing(to_file)
        open(to_file, "w") do io
            write(io, output)
        end
        @info "Results written to $to_file"
    end
    return nothing
end

"""
    mm_to_px(mm, dpi=FIGURE_DPI)

Convert mm to pixels at given DPI.
"""
mm_to_px(mm, dpi=FIGURE_DPI) = round(Int, mm * (1/25.4) * dpi)

"""
    kde_curve(samples; npoints=400)

Evaluate a kernel density estimate of `samples` over an `npoints`-point grid,
returning `(xs, ys)` for plotting.
"""
function kde_curve(samples; npoints=400, boundary=nothing)
    k  = isnothing(boundary) ? kde(samples) : kde(samples, boundary=boundary)
    xs = range(minimum(k.x), maximum(k.x), length=npoints)
    ys = pdf(k, collect(xs))
    return collect(xs), ys
end

"""
    normal_curve(μ, σ; npoints=400, nsigma=4)

Evaluate a `Normal(μ, σ)` pdf over an `npoints`-point grid spanning
`±nsigma` standard deviations, returning `(xs, ys)` for plotting.
"""
function normal_curve(μ, σ; npoints=400, nsigma=4)
    xs = range(μ - nsigma*σ, μ + nsigma*σ, length=npoints)
    ys = pdf.(Normal(μ, σ), xs)
    return collect(xs), ys
end

"""
    ess_rhat_summary(mat)

Effective sample size (bulk/tail) and split-R̂ for a chain-shaped
(iteration × chain) array of values. Works for sampled variables and for any
deterministic transform of them (e.g. derived ages), as long as the chain
structure (not yet flattened across chains) is preserved.
"""
function ess_rhat_summary(mat)
    bulk = ess_rhat(mat; kind=:bulk)
    tail = ess_rhat(mat; kind=:tail)
    return (ess_bulk=bulk.ess, ess_tail=tail.ess, rhat=max(bulk.rhat, tail.rhat))
end

"""
    diagnostics_fig(chain; min_labels, obs_ages, obs_σ, post_vec, post_mat,
                     param_names, param_labels, n_chains)

3-row diagnostic figure:
  Row 1 — trace plot per sampled parameter (`param_names`), with divergent
           draws overlaid in red and R̂ / divergence fraction as in-panel text.
  Row 2 — lower-triangular pairs plot across all sampled parameters: scatter
           of every draw for each pair (divergences highlighted in red), with
           each parameter's own marginal density on the diagonal. Reveals
           joint geometry (funnels, ridges, multimodality) that 1-D trace and
           marginal-density panels cannot show.
  Row 3 — measurement vs. posterior density per physical age (`labels`), ESS
           (bulk/tail) annotated as in-panel text. All panels share a fixed
           x-range (50–250 Ma) and a common y-scale, so relative precision
           across phases is directly comparable.

`post_mat` must hold the same quantities as `post_vec`, but as chain-shaped
(iteration × chain) arrays rather than flattened vectors, so ESS/R̂ can be
computed correctly (flattening loses the per-chain structure needed for that).
"""
function diagnostics_fig(chain; min_labels, obs_ages, obs_σ, post_vec, post_mat,
                          param_names, param_labels, n_chains)
    col_obs      = (:steelblue, 0.8)
    col_post     = (:firebrick, 0.85)
    chain_colors = [:steelblue, :firebrick, :seagreen, :darkorange]
    col_div      = :red

    divergent      = Array(chain[:numerical_error])   # (iteration × chain) Bool
    frac_divergent = sum(divergent) / length(divergent)
    n_params       = length(param_names)

    fig = Figure(size=(mm_to_px(2*190), mm_to_px(2*270)))

    # ── Row 1: Trace plots — divergences overlaid, R̂ + divergence rate ──
    for (col, (pname, plabel)) in enumerate(zip(param_names, param_labels))
        rhat = ess_rhat_summary(Array(chain[pname])).rhat

        ax = Axis(fig[1, col],
                  title  = "Trace: $plabel",
                  xlabel = "Iteration",
                  ylabel = col == 1 ? "Value (Ma)" : "")

        for ch in 1:n_chains
            vals = vec(chain[pname].data[:, ch])
            lines!(ax, 1:length(vals), vals;
                   color=(chain_colors[ch], 0.7), linewidth=0.8)
        end

        text!(ax, 0.98, 0.98;
              text = "R̂=$(round(rhat, digits=3))\ndivergent=$(round(100 * frac_divergent, digits=2))%",
              space = :relative, align = (:right, :top), fontsize = 11)
    end

    # ── Row 2: Pairs plot — joint structure + divergence locations ──────
    Label(fig[2, 1:n_params], "Pairs plot — joint structure of sampled parameters (divergences in red)",
          fontsize=13, font=:bold)

    param_vecs = [vec(chain[pname]) for pname in param_names]
    div_vec    = vec(divergent)

    for i in 1:n_params, j in 1:i
        ax = Axis(fig[2 + i, j],
                  xlabel = i == n_params ? param_labels[j] : "",
                  ylabel = j == 1 ? param_labels[i] : "")

        if i == j
            xs, ys = kde_curve(param_vecs[i])
            lines!(ax, xs, ys; color=:steelblue, linewidth=2)
            band!(ax, xs, zeros(length(xs)), ys; color=(:steelblue, 0.2))
        else
            scatter!(ax, param_vecs[j], param_vecs[i];
                     color=(:steelblue, 0.15), markersize=3, strokewidth=0)
            if any(div_vec)
                scatter!(ax, param_vecs[j][div_vec], param_vecs[i][div_vec];
                         color=col_div, markersize=6)
            end
        end
    end

    # ── Row 3: Measurement vs. posterior, ESS as in-panel text ──────────
    meas_row = 3 + n_params
    axes_age = []
    for (col, (label, μ, σ, samps, mat)) in enumerate(zip(min_labels, obs_ages, obs_σ, post_vec, post_mat))
        s = ess_rhat_summary(mat)

        ax = Axis(fig[meas_row, col],
                  title  = label,
                  xlabel = "Age (Ma)",
                  ylabel = col == 1 ? "Density" : "")

        obs_xs, obs_ys = normal_curve(μ, σ)
        lines!(ax, obs_xs, obs_ys; color=col_obs, linewidth=2.5, linestyle=:dash,
               label="Measured")
        band!(ax, obs_xs, zeros(length(obs_xs)), obs_ys; color=(col_obs[1], 0.2))

        post_xs, post_ys = kde_curve(samps)
        lines!(ax, post_xs, post_ys; color=col_post, linewidth=2.5, label="Posterior")
        band!(ax, post_xs, zeros(length(post_xs)), post_ys; color=(col_post[1], 0.2))

        text!(ax, 0.02, 0.98;
              text = "ESS bulk=$(round(Int, s.ess_bulk))\nESS tail=$(round(Int, s.ess_tail))",
              space = :relative, align = (:left, :top), fontsize = 11)

        xlims!(ax, 100, 180)
        push!(axes_age, ax)
    end
    linkyaxes!(axes_age...)
    axislegend(axes_age[1]; position=:rt, framevisible=false, labelsize=11)

    Label(fig[0, :], "Bayesian Petrochronology — Diagnostics", fontsize=16, font=:bold)

    return fig
end

"""
    meas_v_model_age_fig(; min_labels, colors, obs_ages, obs_σ, post_vec)

Two-panel figure comparing measured age PDFs (left) against posterior age
PDFs (right), one curve per domain, with shared x-axis.
"""
function meas_v_model_age_fig(; min_labels, colors, obs_ages, obs_σ, post_vec)
    fig_compare = Figure(size=(1100, 500))

    ax_meas = Axis(fig_compare[1, 1],
                   title  = "Measured ages",
                   xlabel = "Age (Ma)",
                   ylabel = "Density")

    for (label, μ, σ, color) in zip(min_labels, obs_ages, obs_σ, colors)
        xs, ys = normal_curve(μ, σ)
        lines!(ax_meas, xs, ys; color=color, linewidth=2.5, label="$label ($(μ) ± $(σ) Ma)")
        band!(ax_meas, xs, zeros(length(xs)), ys; color=(color, 0.2))
        vlines!(ax_meas, [μ]; color=color, linestyle=:dash, linewidth=1.2)
    end
    axislegend(ax_meas; position=:lt, framevisible=false, labelsize=11)

    ax_post = Axis(fig_compare[1, 2],
                   title  = "Posterior ages",
                   xlabel = "Age (Ma)",
                   ylabel = "")

    for (label, samps, color) in zip(min_labels, post_vec, colors)
        xs, ys = kde_curve(samps)
        lines!(ax_post, xs, ys; color=color, linewidth=2.5,
               label="$label (median=$(round(median(samps), digits=1)) Ma)")
        band!(ax_post, xs, zeros(length(xs)), ys; color=(color, 0.2))
        vlines!(ax_post, [median(samps)]; color=color, linestyle=:dash, linewidth=1.2)
    end
    axislegend(ax_post; position=:lt, framevisible=false, labelsize=11)

    linkxaxes!(ax_meas, ax_post)

    Label(fig_compare[0, :], "Measured vs. Posterior Age Distributions", fontsize=16, font=:bold)

    return fig_compare
end

"""
    add_ci_bars!(ax, samples_vec, colors; y_pos, central=:median,
                  ci=(0.025, 0.975), dodge_gap=0.02, linewidth=3,
                  markersize=12, whiskerwidth=10)

Overlay horizontal CI bars on `ax` for each vector in `samples_vec`, positioned
around `y_pos` with vertical dodging so bars for different variables don't
overlap. A diamond marker shows the central estimate on each bar.

`y_pos` controls placement — set above or below the density curves as needed.
`central` is `:median` or `:mode` (mode is found as the KDE peak).
`ci` is the lower/upper quantile pair, defaulting to 95%.
"""
function add_ci_bars!(ax, samples_vec, colors; y_pos,
                       central=:median, ci=(0.025, 0.975),
                       dodge_gap=0.05, linewidth=3,
                       markersize=12)
    n       = length(samples_vec)
    offsets = n > 1 ? collect(range(-(n-1)/2, (n-1)/2, length=n)) .* dodge_gap : [0.0]

    for (samps, color, offset) in zip(samples_vec, colors, offsets)
        y  = y_pos + offset
        lo = quantile(samps, ci[1])
        hi = quantile(samps, ci[2])

        cval = if central == :median
            median(samps)
        else
            xs, ys = kde_curve(samps)
            xs[argmax(ys)]
        end

        rangebars!(ax, [y], [lo], [hi];
                   direction=:x, color=color,
                   linewidth=linewidth)
        scatter!(ax, [cval], [y];
                 color=color, marker=:diamond, markersize=markersize)
    end
end
