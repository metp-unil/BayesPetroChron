μ_ages = (
    grt_core = 136.7,
    grt_rim  = 139.7,
    mnz      = 138.8,
    zrn      = 137.2,
    rt       = 121.5,
)

σ_ages = (
    grt_core = 20.1,
    grt_rim  = 8.7,
    mnz      = 3.4,
    zrn      = 2.7,
    rt       = 2.5,
)

obs_ages = [μ_ages.grt_core, μ_ages.grt_rim, μ_ages.mnz, μ_ages.zrn, μ_ages.rt]
obs_σ    = [σ_ages.grt_core, σ_ages.grt_rim, σ_ages.mnz, σ_ages.zrn, σ_ages.rt] .* 0.5    # reported as 2σ

min_labels = ["Garnet core", "Garnet rim", "Monazite", "Zircon", "Rutile"]
# alt colours
min_colors = ["#79132E", "#DB5C7E", "#2C895E", "#83892C", "#003D6E"]
# Gutos colours
min_colors = ["#5201A3", "#FADA24", "#D97A3A", "#007C83", "#3B5B92"]