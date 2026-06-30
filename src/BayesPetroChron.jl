module BayesPetroChron

using CairoMakie
using KernelDensity
using Distributions
using Statistics
using MCMCDiagnosticTools

export diagnostics_fig, meas_v_model_age_fig

include("plotting.jl")

end
