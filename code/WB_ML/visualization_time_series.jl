using Trixi
using DelimitedFiles
using DataFrames
using CairoMakie
using LaTeXStrings

####################################################################################################
# Plotting routine to visualize the time series of the lake-at-rest error and the entropy 
# time-derivative for a lake-at-rest setup with and without perturbation
####################################################################################################

####################################################################################################
# Run elixirs with and without perturbation

trixi_include("elixir_shallowwater_multilayer_well_balanced_wet_dry.jl",
              perturbation = true,
              output_directory = "wb_perturbation")
trixi_include("elixir_shallowwater_multilayer_well_balanced_wet_dry.jl",
              perturbation = false,
              output_directory = "wb")

####################################################################################################
# Read data from analysis.dat into Dataframes

mat, head = readdlm(joinpath("wb", "analysis.dat"), header = true)
df_wb = DataFrame(mat, vec(head))

mat, head = readdlm(joinpath("wb_perturbation", "analysis.dat"), header = true)
df_wb_perturbation = DataFrame(mat, vec(head))

####################################################################################################
# Create plot

with_theme(theme_latexfonts()) do
    # Setup figure
    f = Figure(size = (550, 550 / 2.5))
    ax_left = Axis(f[1, 1], xlabel = "time", ylabel = L"\bar{S}_t", yscale = Makie.log10,
                   yreversed = true,
                   yticks = ([1e-30, 1e-20, 1e-10, 1e-0],
                             [L"-10^{-30}", L"-10^{-20}", L"-10^{-10}", L"-10^{0}"]))
    ax_right = Axis(f[1, 2], xlabel = "time",
                    ylabel = L"\frac{1}{|\Omega|}\int|H_1-H_1(0)|d\Omega",
                    yscale = Makie.log10,
                    yticks = ([1e-16, 1e-12, 1e-8, 1e-4, 1e-0],
                              [L"10^{-16}", L"10^{-12}", L"10^{-8}", L"10^{4}", L"10^{0}"]))
    Makie.ylims!(ax_right, 1e-16, 1e-0)
    Makie.ylims!(ax_left, 1e-35, 1e-0)

    # Plot lake-at-rest error and entropy time-derivative
    lines!(ax_left, df_wb.time[2:3:end], abs.(df_wb.dsdu_ut)[2:3:end],
           label = "steady-state")
    lines!(ax_left, df_wb_perturbation.time[2:3:end],
           abs.(df_wb_perturbation.dsdu_ut)[2:3:end], label = "perturbed")
    lines!(ax_right, df_wb.time[2:3:end], df_wb.var"|H0-(h+b)|"[2:3:end],
           label = "steady-state")
    lines!(ax_right, df_wb_perturbation.time[2:3:end],
           df_wb_perturbation.var"|H0-(h+b)|"[2:3:end], label = "perturbed")

    # Add legends
    Legend(f[2, 1:2], ax_left, labelsize = 11, orientation = :horizontal,
           framevisible = :false)

    # Adjust row and column gaps
    rowgap!(f.layout, 0)
    colgap!(f.layout, 50)

    save("time_series_wb_ml.pdf", f)
end
