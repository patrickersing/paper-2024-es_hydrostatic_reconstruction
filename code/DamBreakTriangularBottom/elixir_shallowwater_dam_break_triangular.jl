using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the shallow water equations

equations = ShallowWaterMultiLayerEquations1D(gravity_constant = 9.812, rhos = 1.0,
                                              threshold_desingularization = 1e-8)

"""
    initial_condition_dam_break_triangular(x, t, equations:: ShallowWaterMultiLayerEquations1D)

Initial condition for an experimental test case of a dam break over triangular bottom topography
with wetting and drying. 
For this test case experimental data for the water height is available at several gauge locations.

More information about this test can be found in the papers:
- S. Gu et al. (2017)
  SWE-SPHysics Simulation of Dam Break Flows at South-Gate Gorges Reservoir
  [DOI: 10.3390/w9060387](https://doi.org/10.3390/w9060387)
- J.G. Zhou, D.M. Causon, C.G. Mingham and D.M. Ingram (2004)
  Numerical Prediction of Dam-Break Flows in General Geometries with Complex Bed Topography
  [DOI: 10.1061/(ASCE)0733-9429(2004)130:4(332)](https://doi.org/10.1061/(ASCE)0733-9429(2004)130:4(332))
"""
function initial_condition_dam_break_triangular(x, t,
                                                equations::ShallowWaterMultiLayerEquations1D)
    # Set background values                                             
    b = 0.0
    h = 0.0
    v = 0.0

    if x[1] <= 15.5
        h = 0.75
    elseif 25.5 < x[1] && x[1] <= 28.5
        b = (x[1] - 25.5) * 0.4 / 3.0
    elseif x[1] > 28.5 && x[1] < 31.5
        b = 0.4 - (x[1] - 28.5) * 0.4 / 3.0
    end

    H = h + b

    if x[1] > 28.5
        H = max(H, 0.15)
    end

    # It is mandatory to shift the water level at dry areas to make sure the water height h
    # stays positive. The system would not be stable for h set to a hard 0 due to division by h in
    # the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    # with a default value of 5*eps() ≈ 1e-13 in double precision, is set in the constructor above
    # for the ShallowWaterEquations and added to the initial condition if h = 0.
    # This default value can be changed within the constructor call depending on the simulation setup.
    H = max(H, b + equations.threshold_limiter)
    return prim2cons(SVector(H, v, b), equations)
end

# Source term for manning friction
@inline function source_term_manning_friction(u, x, t,
                                              equations::ShallowWaterMultiLayerEquations1D)
    h, hv, _ = u

    # Friction Coefficient
    n = 0.0125
    # Desingularization
    h = (2.0 * h) / (h^2 + max(h^2, 1e-8))

    return SVector(0.0, -equations.gravity * n^2 * h^(7 / 3) * abs(hv) * hv, 0.0)
end

initial_condition = initial_condition_dam_break_triangular

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                              hydrostatic_reconstruction_ersing_etal))

basis = LobattoLegendreBasis(4)

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = waterheight_pressure)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Create the TreeMesh for the domain [0, 38]

coordinates_min = 0.0
coordinates_max = 38.0

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 7,
                n_cells_max = 20000,
                periodicity = false)

boundary_condition = boundary_condition_slip_wall

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition,
                                    source_terms = source_term_manning_friction)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 40.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 5000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     extra_analysis_errors = (:conservation_error,),
                                     save_analysis = false)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 1000,
                                     save_initial_solution = true,
                                     save_final_solution = true)

# Save time_series data at the gauge locations                             
time_series = TimeSeriesCallback(semi, [(19.5), (25.5), (28.5), (35.5)];
                                 interval = 1,
                                 solution_variables = cons2cons,
                                 filename = "tseries.h5")

stepsize_callback = StepsizeCallback(cfl = 0.7)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        time_series,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(stage_limiter!); dt = 1.0,
            ode_default_options()..., callback = callbacks, adaptive = false);

summary_callback() # print the timer summary
