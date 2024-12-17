using Trixi
using TrixiShallowWater
using Printf: @printf, @sprintf
using CairoMakie
using LaTeXStrings

####################################################################################################
# Plot routine to visualize the approximated solution for the 2D Parabolic bowl test case at 
# different times.
####################################################################################################
T = 2 * pi / sqrt(2*9.81*0.1) # Period of the oscillation

# Compute numerical approximation at t=4T/3
trixi_include("elixir_shallowwater_multilayer_parabolic_bowl.jl", tspan = (0.0, 4T/3))
pd = PlotData1D(sol, slice = :x, point = (0,0))

# Exact solution
x_exact = LinRange(-2, 2, 500)
H_exact = [initial_condition_parabolic_bowl((x, 0), tspan[2], equations)[1] + initial_condition_parabolic_bowl((x,0), tspan[2], equations)[4] for x in x_exact]
v_exact = [cons2prim(initial_condition_parabolic_bowl((x, 0), tspan[2], equations), equations)[2] for x in x_exact]
w_exact = [cons2prim(initial_condition_parabolic_bowl((x, 0), tspan[2], equations), equations)[3] for x in x_exact]
hv_exact = [initial_condition_parabolic_bowl((x, 0), tspan[2], equations)[2] for x in x_exact]
hw_exact = [initial_condition_parabolic_bowl((x, 0), tspan[2], equations)[3] for x in x_exact]

# Plot water height and bottom topography
with_theme(theme_latexfonts()) do
  f = Figure(size = (400, 300))
  ax = Axis(f[1,1], xlabel = "x", ylabel = "height", limits = (nothing, (nothing, 0.5)))

  lines!(ax, pd.x .- 2, pd.data[:,1]; linewidth = 2.0, linestyle = :dash, color = :black, markersize=6, label = "numerical")
  lines!(ax, x_exact, H_exact; linewidth = 1.0, linestyle = :solid, color = :black, label = "exact")
  lines!(ax, pd.x .- 2, pd.data[:,4]; linewidth = 1.0, linestyle = :dot, color =:grey, label = "bottom")

  # Add legend
  axislegend(ax, orientation = :horizontal, framevisible = :true, location = :best)
  display(f)

  save("visualization_parabolic_bowl_2d_N3_1T3.pdf", f)
end

# Momentum plots
with_theme(theme_latexfonts()) do
  f = Figure(size = (400, 300))
  ax = Axis(f[1,1], xlabel = "x", ylabel = "momentum", limits = (nothing, (nothing, 0.1)))

  lines!(ax, pd.x .- 2, pd.data[:,2].*(pd.data[:,1].-pd.data[:,4]); linewidth = 2.0, linestyle = :dash, color=Cycled(1), label = L"hv_{num}")
  lines!(ax, x_exact, hv_exact; linewidth = 1.0, linestyle = :solid, color=Cycled(1), label = L"hv_{ref}")
  lines!(ax, pd.x .- 2, pd.data[:,3].*(pd.data[:,1].-pd.data[:,4]); linewidth = 2.0, linestyle = :dash, color=Cycled(2), label = L"hw_{num}")
  lines!(ax, x_exact, hw_exact; linewidth = 1.0, linestyle = :solid, color=Cycled(2), label = L"hw_{ref}")

  # Add legend
  axislegend(ax, orientation = :horizontal, framevisible = :true, location = :best)
  display(f)
  save("visualization_parabolic_bowl_2d_N3_1T3_mom.pdf", f)
end

# Compute numerical approximation at t=2t
trixi_include("elixir_shallowwater_multilayer_parabolic_bowl.jl", tspan = (0.0, 2T))
pd = PlotData1D(sol, slice = :x, point = (0,0))

# Exact solution
x_exact = LinRange(-2, 2, 500)
H_exact = [initial_condition_parabolic_bowl((x, 0), tspan[2], equations)[1] + initial_condition_parabolic_bowl((x,0), tspan[2], equations)[4] for x in x_exact]
v_exact = [cons2prim(initial_condition_parabolic_bowl((x, 0), tspan[2], equations), equations)[2] for x in x_exact]
w_exact = [cons2prim(initial_condition_parabolic_bowl((x, 0), tspan[2], equations), equations)[3] for x in x_exact]
hv_exact = [initial_condition_parabolic_bowl((x, 0), tspan[2], equations)[2] for x in x_exact]
hw_exact = [initial_condition_parabolic_bowl((x, 0), tspan[2], equations)[3] for x in x_exact]

# Plot water height and bottom topography
with_theme(theme_latexfonts()) do
  f = Figure(size = (400, 300))
  ax = Axis(f[1,1], xlabel = "x", ylabel = "height", limits = (nothing, (nothing, 0.5)))

  lines!(ax, pd.x .- 2, pd.data[:,1]; linewidth = 2.0, linestyle = :dash, color = :black, markersize=6, label = "numerical")
  lines!(ax, x_exact, H_exact; linewidth = 1.0, linestyle = :solid, color = :black, label = "exact")
  lines!(ax, pd.x .- 2, pd.data[:,4]; linewidth = 1.0, linestyle = :dot, color =:grey, label = "bottom")

  # Add legend
  axislegend(ax, orientation = :horizontal, framevisible = :true, location = :best)
  display(f)

  save("visualization_parabolic_bowl_2d_N3_2T.pdf", f)
end

# Momentum plots
with_theme(theme_latexfonts()) do
  f = Figure(size = (400, 300))
  ax = Axis(f[1,1], xlabel = "x", ylabel = "momentum", limits = (nothing, (nothing, 0.1)))

  lines!(ax, pd.x .- 2, pd.data[:,2].*(pd.data[:,1].-pd.data[:,4]); linewidth = 2.0, linestyle = :dash, color=Cycled(1), label = L"hv_{num}")
  lines!(ax, x_exact, hv_exact; linewidth = 1.0, linestyle = :solid, color=Cycled(1), label = L"hv_{ref}")
  lines!(ax, pd.x .- 2, pd.data[:,3].*(pd.data[:,1].-pd.data[:,4]); linewidth = 2.0, linestyle = :dash, color=Cycled(2), label = L"hw_{num}")
  lines!(ax, x_exact, hw_exact; linewidth = 1.0, linestyle = :solid, color=Cycled(2), label = L"hw_{ref}")

  # Add legend
  axislegend(ax, orientation = :horizontal, framevisible = :true, location = :best)
  display(f)
  save("visualization_parabolic_bowl_2d_N3_2T_mom.pdf", f)
end