using Trixi
using TrixiShallowWater
using Printf: @printf, @sprintf
using CairoMakie
using LaTeXStrings

####################################################################################################
# Plot routine to visualize the experiment setup with gauge locations, initial condition and 
# temporal evolution of the water height.
####################################################################################################

# Plot initial condition and bottom topography
####################################################################################################
#= Due to the use of LGL nodes, the projected initial condition will not have a true discontinuity.
Therefore, to visualize the experimental setup we only extract the node positions and then 
directly plot the initial condition instead of the projection. =#

# Run for one time step to obtain the node positions
trixi_include("elixir_shallowwater_dam_break_triangular.jl", tspan = (0.0, 1e-10))
pd = PlotData1D(sol)

# Compute the true initial condition at the node positions
init_data = initial_condition_dam_break_triangular.(pd.x, 0.0, equations)
h0 = [u[1] for u in init_data]
b = [u[3] for u in init_data]
H0 = h0 + b

with_theme(theme_latexfonts()) do
    f = Figure(size = (550, 550 / 2.5))
    # Creat additional labels to indicate gauge locations
    ax = Axis(f[1, 1], xlabel = "x / m", ylabel = "waterheight / m",
              limits = (0, 38, 0.0, 1.2),
              xticks = ([0, 15.5, 15.5 + 4.0, 15.5 + 10.0, 15.5 + 13.0, 15.5 + 20.0, 38],
                        ["0.0", "15.5", "G4", "G10", "G13", "G20", "38.0"]))

    # Plot temporal evolution of the waterheight at different times
    t = [5, 10, 20, 30]
    l = []  # Legend entries need to be reordered
    for (i, t_i) in enumerate(t)
        trixi_include("elixir_shallowwater_dam_break_triangular.jl", tspan = (0.0, t_i))
        pd = PlotData1D(sol)
        push!(l,
              lines!(ax, pd.x, pd.data[:, 1], label = @sprintf("t = %.1fs", t_i), color = i,
                     colormap = :tab20c, colorrange = (1, 20), linewidth = 2.0))
    end

    # Plot initial condition / bottom
    l0 = lines!(ax, pd.x, H0, linewidth = 2.0, linestyle = :dash, color = :black,
                label = @sprintf("t = %.1fs", 0.0))
    band!(ax, pd.x, 0.0, b, color = :gray95)  # Set color for bottom topography
    lines!(ax, pd.x, b, linewidth = 1.5, color = :black, linestyle = :solid)

    # Add legend
    axislegend(ax, [l0, l...], [L"t=0.0s", L"t=5.0s", L"t=10.0s", L"t=20.0s", L"t=30.0s"],
               labelsize = 11, orientation = :horizontal)

    save("visualization_dam_break_triangular.pdf", f)
end
