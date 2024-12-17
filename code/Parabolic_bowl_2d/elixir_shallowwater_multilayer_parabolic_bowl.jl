
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations for a single layer

equations = ShallowWaterMultiLayerEquations2D(gravity_constant = 9.81, rhos = 1.0,
                                              threshold_desingularization = 1e-8)

"""
    initial_condition_parabolic_bowl(x, t, equations:: ShallowWaterMultiLayerEquations2D)

Well-known initial condition to test the [`hydrostatic_reconstruction_ersing_etal`](@ref) and its
wet-dry mechanics. This test has an analytical solution. The initial condition is defined by the
analytical solution at time t=0. The bottom topography defines a bowl and the water level is given
by an oscillating lake.

The original test and its analytical solution were first presented in
- William C. Thacker (1981)
  Some exact solutions to the nonlinear shallow-water wave equations
  [DOI: 10.1017/S0022112081001882](https://doi.org/10.1017/S0022112081001882).

The particular setup below is taken from Section 6.2 of
- Niklas Wintermeyer, Andrew R. Winters, Gregor J. Gassner and Timothy Warburton (2018)
  An entropy stable discontinuous Galerkin method for the shallow water equations on
  curvilinear meshes with wet/dry fronts accelerated by GPUs
  [DOI: 10.1016/j.jcp.2018.08.038](https://doi.org/10.1016/j.jcp.2018.08.038).
"""
function initial_condition_parabolic_bowl(x, t, equations::ShallowWaterMultiLayerEquations2D)
    a = 1.0
    h_0 = 0.1
    sigma = 0.5
    ω = sqrt(2 * equations.gravity * h_0) / a

    b = h_0 * ((x[1])^2 + (x[2])^2) / a^2

    H = sigma * h_0 / a^2 * (2 * x[1] * cos(ω * t) + 2 * x[2] * sin(ω * t) - sigma) + h_0

    if (H-b) > 0.0
    v1 = -sigma * ω * sin(ω * t)
    v2 = sigma * ω * cos(ω * t)
    else
      v1 = 0.0
      v2 = 0.0
    end
    # It is mandatory to shift the water level at dry areas to make sure the water height h
    # stays positive. The system would not be stable for h set to a hard 0 due to division by h in
    # the computation of velocity, e.g., (h v1) / h. Therefore, a small dry state threshold
    # with a default value of 500*eps() ≈ 1e-13 in double precision, is set in the constructor above
    # for the ShallowWaterEquationsWetDry and added to the initial condition if h = 0.
    # This default value can be changed within the constructor call depending on the simulation setup.
    H = max(H, b + equations.threshold_limiter)
    return prim2cons(SVector(H, v1, v2, b), equations)
end

initial_condition = initial_condition_parabolic_bowl
###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                              hydrostatic_reconstruction_ersing_etal))

basis = LobattoLegendreBasis(3)

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
# Create the StructuredMesh for the domain [-2, 2]^2

coordinates_min = (-2.0, -2.0)
coordinates_max = (2.0, 2.0)

cells_per_dimension = (100, 100)

mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solvers, callbacks etc.

T = 2 * pi / sqrt(2*9.81*0.1) # period of the oscillation

tspan = (0.0, 2*T)

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false,
                                     extra_analysis_errors = (:conservation_error,))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = T/4,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.7)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution, stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(stage_limiter!); dt = 1.0,
            ode_default_options()..., callback = callbacks, adaptive = false);

summary_callback() # print the timer summary