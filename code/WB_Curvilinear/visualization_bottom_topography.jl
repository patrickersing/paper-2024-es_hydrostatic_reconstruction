using Trixi
using DataFrames
using Statistics
using DelimitedFiles
using Printf
using CairoMakie
using LaTeXStrings

####################################################################################################
# Ploting routine create a contour plot of the bottom topography and visualize the curvilinear mesh
####################################################################################################

# Run elixir
trixi_include("elixir_shallowwater_multilayer_well_balanced_wet_dry.jl")

# Create plots
with_theme(theme_latexfonts()) do
    f = Figure(size = (300, 200))

    # Create contour plot of bottom topoography
    layout_left = f[1, 1]
    ax_left = Axis(layout_left[1, 1], aspect = 1, xlabel = "x", ylabel = "y",
                   xticks = ([0.0, sqrt(2) / 2, sqrt(2)],
                             [L"0", L"\sqrt{2}/2", L"\sqrt{2}"]),
                   yticks = ([0.0, sqrt(2) / 2, sqrt(2)],
                             [L"0", L"\sqrt{2}/2", L"\sqrt{2}"]))
    pd = PlotData2D(sol)
    cont_plot = Makie.plot!(ax_left, pd["b"], colormap = (:viridis, 1.0),
                            show_colorbar = false, limits = (0.0, sqrt(2));
                            plot_mesh = true)
    Makie.Colorbar(layout_left[1, 2], cont_plot, label = "b", vertical = true)
    tightlimits!(ax_left)

    save("bottom_topography_wb_curvilinear.pdf", f)
end
