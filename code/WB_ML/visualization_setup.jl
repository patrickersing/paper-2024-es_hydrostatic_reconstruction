using Trixi
using TrixiShallowWater
using Printf: @printf, @sprintf
using LaTeXStrings
using CairoMakie

####################################################################################################
# Plotting routine to visualize the perturbed initial condition together with a temporal evolution
# of the layer heights
####################################################################################################

with_theme(theme_latexfonts()) do
    # Setup figure and axis
    f = Figure(size = (550, 550 / 2.5))
    ax = Axis(f[1, 1], xlabel = "x / m", ylabel = "waterheight / m",
              limits = (0, 100, 0.0, 3.5))

    # Calculate elixir up to different times to plot the temporal evolution
    t = [1e-16, 10, 100, 1000] # times to plot
    ls = [:solid, :dash, :dot, :dashdot] # linestyles

    for (i, t_i) in enumerate(t)
        # run elixir to obtain solution data
        trixi_include("elixir_shallowwater_multilayer_well_balanced_wet_dry.jl",
                      tspan = (0.0, t_i))
        pd = PlotData1D(sol)

        # Plot layer heights
        lines!(ax, pd.x, pd.data[:, 1], label = @sprintf("t = %.1fs", t_i), color = i,
               colormap = (:tab20c), linewidth = 1.5, colorrange = (1, 20),
               linestyle = (ls[i], :dense))
        lines!(ax, pd.x, pd.data[:, 2], color = i, colormap = (:tab20c),
               colorrange = (1, 20), linewidth = 1.5, linestyle = (ls[i], :dense))
    end

    # Plot bottom topography
    trixi_include("elixir_shallowwater_multilayer_well_balanced_wet_dry.jl",
                  tspan = (0.0, 1e-16))
    pd = PlotData1D(sol)
    lines!(ax, pd.x, pd.data[:, end], linewidth = 2.0, color = :black)
    band!(ax, pd.x, 0.0, pd.data[:, end], color = :gray95)  # Set color for bottom topography

    # Add legend
    axislegend(ax, labelsize = 11, orientation = :horizontal)

    save("visualization_setup_wb_ml.pdf", f)
end
