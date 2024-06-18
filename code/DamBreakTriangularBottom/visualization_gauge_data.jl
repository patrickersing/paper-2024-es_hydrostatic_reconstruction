
using Trixi
using TrixiShallowWater
using CSV
using DataFrames
using CairoMakie

####################################################################################################
# Create time series plot for different gauge locations
####################################################################################################

with_theme(theme_latexfonts()) do
    # Define plot positions and gauge names
    pos = ((1, 1), (1, 2), (2, 1), (2, 2))
    names = ["G4", "G10", "G13", "G20"]

    # Create layout
    f = Figure(size = (550, 400))
    ax_n = Vector{Any}(undef, 4)
    # Workaround to avoid multiple labels
    for i in 1:4
        if i == 1 || i == 3
            ax_n[i] = Axis(f[pos[i]...], ylabel = "waterheight / m", title = names[i],
                           limits = (0, 40, 0.0, 0.8))
        else
            ax_n[i] = Axis(f[pos[i]...], title = names[i], limits = (0, 40, 0.0, 0.8))
        end
    end
    ax_n[3].xlabel = "time / s"
    ax_n[4].xlabel = "time / s"

    # Adjust row and column gaps
    colgap!(f.layout, 50)
    rowgap!(f.layout, 10)

    linkaxes!(ax_n...)

    # Run elixir to obtain solution data
    include("elixir_shallowwater_dam_break_triangular.jl")
    time_series_new = time_series

    # Loop over gauge locations and plot data
    for (i, name) in enumerate(names)
        # Plot experimental data
        pd_exp = CSV.read(joinpath("Reference", name * "_Experimental.csv"), DataFrame)
        Makie.scatter!(ax_n[i], pd_exp.X, pd_exp.Y, label = "Experiment", color = 2,
                       colormap = (:tab10, 0.5), colorrange = (1, 10), markersize = 7,
                       strokewidth = 1, strokecolor = (:black, 0.5))

        # Plot time_series
        pd = PlotData1D(time_series_new, i)
        lines!(ax_n[i], pd.x, pd.data[:, 1], label = "New HR", linewidth = 1.5)

        # Add legend
        Legend(f[3, 1:2], ax_n[1], orientation = :horizontal, framevisible = :false)
    end

    save("time_series_dam_break_triangular.pdf", f)
end
