using Trixi
using TrixiShallowWater
using Printf
using DelimitedFiles
using LaTeXStrings
using DataFrames
using CairoMakie

# Run test case
####################################################################################################
# set interval for polydeg
n_start = 6
n_end = 24

# Loop over polydeg and save analysis files
for i in n_start:n_end
    trixi_include("elixir_shallowwater_multilayer_convergence.jl", polydeg = i,
                  output_directory = joinpath("ES", @sprintf("P%i", i)))
end

# Organize data and create plot
####################################################################################################
# Create Vectors to save plotting data
L2_h1 = Vector{Float64}()
L2_h2 = Vector{Float64}()
L2_h3 = Vector{Float64}()
N = Vector{Int16}()
dt = 0

cd("ES")
dir_names = readdir()
# Iterate over folders for different polydeg
for name in dir_names
    m = match(r"P(\d+)", name)
    # If Folder "PX" found read analysis.dat file
    if (m !== nothing)
        # Read from analysis.dat file
        analysis_file = joinpath(name, "analysis.dat")
        mat, head = readdlm(analysis_file, header = true)
        df = DataFrame(mat, vec(head))

        # Save L2_Errors and polydeg N
        push!(L2_h1, df.l2_h1[end])
        push!(L2_h2, df.l2_h2[end])
        push!(L2_h3, df.l2_h3[end])
        push!(N, parse(Int16, m.captures[1]))
        dt = Int(round(1 / df.dt[1]))
    end
end
cd("..")

# Read into Dataframe and sort
df_conv = DataFrame(N_pol = N, L2_h1 = L2_h1, L2_h2 = L2_h2, L2_h3 = L2_h3)
df_conv = sort(df_conv, :N_pol)

# Plot data
####################################################################################################
with_theme(theme_latexfonts()) do
    f = Figure(size = (400, 400 / 1.618))
    ax = Axis(f[1, 1], xlabel = L"N", ylabel = L"L_2\text{ Error}", yscale = Makie.log10,
              yminorticksvisible = true, yminorticks = IntervalsBetween(8))

    Makie.scatterlines!(ax, df_conv.N_pol, df_conv.L2_h1, label = L"h_1",
                        marker = (:circle))
    Makie.scatterlines!(ax, df_conv.N_pol, df_conv.L2_h2, label = L"h_2", marker = (:rect))
    Makie.scatterlines!(ax, df_conv.N_pol, df_conv.L2_h3, label = L"h_3",
                        marker = (:utriangle))

    axislegend(ax, labelsize = 11, orientation = :horizontal)

    save("spectral_convergence_plot.pdf", f)
end
