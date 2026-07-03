μ_ages = (
    grt_core = 136.7,
    grt_rim  = 139.7,
    mnz      = 138.8,
#    zrn      = 137.2,      # Frankfurt date
    zrn      = 136.4,
    rt       = 121.5,
)

σ_ages = (
    grt_core = 20.1,
    grt_rim  = 8.7,
    mnz      = 3.4,
#    zrn      = 1.2,        # Frankfurt date
    zrn      = 1.3,
    rt       = 2.5,
)

obs_ages = [μ_ages.grt_core, μ_ages.grt_rim, μ_ages.mnz, μ_ages.zrn, μ_ages.rt]
obs_σ    = [σ_ages.grt_core, σ_ages.grt_rim, σ_ages.mnz, σ_ages.zrn, σ_ages.rt] .* 0.5    # reported as 2σ

min_labels = ["Garnet core", "Garnet rim", "Monazite", "Zircon", "Rutile"]
# min_colors = [:firebrick, :orange, :seagreen, :steelblue, :purple]
min_colors = ["#643174", "#EDCF3D", "#F48C83", "#CE8D67", "#003D6E"]
min_colors = ["#79132E", "#DB5C7E", "#2C895E", "#83892C", "#003D6E"]
