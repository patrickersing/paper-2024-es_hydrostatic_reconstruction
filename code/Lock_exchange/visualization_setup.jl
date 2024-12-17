using Trixi
using TrixiShallowWater
using Printf: @printf, @sprintf
using CairoMakie
using LaTeXStrings

####################################################################################################
# Plot routine to visualize the final solution for poldeg=1 and polydeg = 2 together with an 
# indicator for the loss of hyperbolicity.
####################################################################################################

# Compute numerical approximation for N=1
trixi_include("elixir_shallowwater_multilayer_lock_exchange.jl", polydeg=1)
pd_N1 = PlotData1D(sol)

# Compute hyperbolicity indicator
ind_hyp_N1 = (pd_N1.data[:,3] - pd_N1.data[:,4]).^2 ./ (9.81*(1 - 0.98)*(pd_N1.data[:,1] - pd_N1.data[:,5]))
H1_N1 = pd_N1.data[:,1]
H2_N1 = pd_N1.data[:,2]

# Compute numerical approximation for N=2
trixi_include("elixir_shallowwater_multilayer_lock_exchange.jl", polydeg=2)
pd_N2 = PlotData1D(sol)

# Compute hyperbolicity indicator
ind_hyp_N2 = (pd_N2.data[:,3] - pd_N2.data[:,4]).^2 ./ (9.81*(1 - 0.98)*(pd_N2.data[:,1] - pd_N2.data[:,5]))
H1_N2 = pd_N2.data[:,1]
H2_N2 = pd_N2.data[:,2]


with_theme(theme_latexfonts()) do
    # Setup the figure and axes
    f = Figure(size = (600, 350))
    ax_bottom = Axis(f[2,1], xlabel = "x", ylabel = "waterheight")
    ax_top = Axis(f[1,1], ylabel = L"I_{hyp}")

    rowsize!(f.layout, 1, Relative(0.4))    # adjust spacing between the two axes

    # Set colors
    H1_color = Cycled(1)
    H2_color = Cycled(2)
    ind_color = :black

    # Plot upper layer height
    lines!(ax_bottom, pd_N2.x, H1_N2, label = L"H_1\;(N = 2)", color = H1_color, linestyle=:solid, linewidth = 1.0)
    lines!(ax_bottom, pd_N1.x, H1_N1, label = L"H_1\;(N = 1)", color = H1_color, linestyle=:dash, linewidth=1.0)

    # Plot lower layer height
    lines!(ax_bottom, pd_N2.x, H2_N2, label = L"H_2\;(N=2)", color = H2_color, linestyle=:solid, linewidth = 1.0)
    lines!(ax_bottom, pd_N1.x, H2_N1, label = L"H_2\;(N=1)", color = H2_color, linestyle=:dash, linewidth=1.0)

    # Plot bottom topography
    lines!(ax_bottom, pd_N1.x, pd_N1.data[:,5], label = L"b", color = :black, linestyle = :solid, linewidth = 1.0)

    # Plot hyperbolicity indicator
    lines!(ax_top, pd_N1.x, ind_hyp_N1, label = L"N=1", color = ind_color, linestyle = :dash, linewidth=1.0)
    lines!(ax_top, pd_N2.x, ind_hyp_N2, label = L"N=2", color = ind_color, linestyle=:solid, linewidth=1.0)

    # Line to denote transition location for hyperbolicity loss
    vlines!(ax_top, -0.21, color = :red, linestyle = :solid, linewidth=1)
    vlines!(ax_bottom, -0.21, color = :red, linestyle = :solid, linewidth=1)

    # Add Legend
    Legend(f[2,2], ax_bottom, labelsize = 11, orientation = :vertical, location = :best)
    Legend(f[1,2], ax_top, labelsize = 11, orientation = :vertical, location = :best)

    # Save figure
    save("visualization_lock_exchange.pdf", f)
end
