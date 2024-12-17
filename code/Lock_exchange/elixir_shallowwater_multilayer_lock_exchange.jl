
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater
using StaticArrays: MVector

###############################################################################
# Semidiscretization of the multilayer shallow water equations for a lock exchange flow with two
# layers and an initial dam break over a Gaussian shaped bottom topography

equations = ShallowWaterMultiLayerEquations1D(gravity_constant = 9.81,
                                              rhos = (0.98, 1.0))


"""
    function initial_condition_lock_exchange(x, t, equations::ShallowWaterMultiLayerEquations1D)

Initial condition for a lock exchange problem with two layers over a Gaussian shaped bottom topography.
The test is initialized with lighter fluid on the left side and heavier fluid on the right side.
Both layers are at rest and separated by a discontinuity in the layer heights. 

This test case was considered in:
- A. Kurganov, G. Petrova (2009)
  Central-upwind schemes for two-layer shallow water equations
  [DOI: 10.1137/08071909](https://doi.org/10.1137/08071909)
- F. Benkhaldoun, S. Sari, M. Seaid (2014)
  A simple multi-layer finite volume solver for density-driven shallow water flows
  [DOI: 10.1016/j.matcom.2013.04.016](https://doi.org/10.1016/j.matcom.2013.04.016)
- E.D. Fernández-Nieto, M.J. Castro Díaz, C. Parés (2011)
  On an intermediate field capturing Riemann solver based on a parabolic viscosity matrix for the two-layer shallow water system
  [DOI: 10.1007/s10915-011-9465-7](https://doi.org/10.1007/s10915-011-9465-7)
"""
function initial_condition_lock_exchange(x, t, equations::ShallowWaterMultiLayerEquations1D)
    b = exp(-x[1]^2) - 2.0

    # Set the discontinuity
    if x[1] <= 0
        H = [0.0, b]
    else
        # Right side of the domain is dry
        H = [0.0, 0.0]
    end

    # It is mandatory to shift the water level at dry areas to make sure the water height h
    # stays positive. The system would not be stable for h set to a hard 0 due to division by h in
    # the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    # with a default value of 5*eps() ≈ 1e-15 in double precision, is set in the constructor above
    # for the ShallowWaterMultiLayerEquations1D and added to the initial condition if h = 0.
    # This default value can be changed within the constructor call depending on the simulation setup.
    for i in reverse(eachlayer(equations))
        if i == nlayers(equations)
            H[i] = max(H[i], b + equations.threshold_limiter)
        else
            H[i] = max(H[i], H[i + 1] + equations.threshold_limiter)
        end
    end

    # Initialize zero velocity
    v = zero(H)

    return prim2cons(SVector(H..., v..., b), equations)
end

function boundary_condition_lock_exchange(u_inner, orientation_or_normal, direction,
    x, t, surface_flux_function,
    equations::ShallowWaterMultiLayerEquations1D)

    # Create the "external" boundary solution state
    h = Trixi.waterheight(u_inner, equations)
    hv = TrixiShallowWater.momentum(u_inner, equations)
    b = u_inner[end]

    # Set the boundary state
    if x[1] > 0
        if hv[1] > eps()
            u_boundary = SVector(h..., hv[1], -hv[1], b)
        else
            u_boundary = SVector(h..., -hv..., b)
        end
    elseif hv[2] < eps()
        u_boundary = SVector(h[1], h[2], -hv[2], hv[2], b)
    else
        u_boundary = SVector(h..., -hv..., b)
    end

    # Calculate the boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        f = surface_flux_function(u_inner, u_boundary, orientation_or_normal, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        f = surface_flux_function(u_boundary, u_inner, orientation_or_normal, equations)
    end

    return f
end

initial_condition = initial_condition_dam_break
boundary_condition = boundary_condition_lock_exchange

###############################################################################
# Get the DG approximation space

volume_flux = (flux_ersing_etal, flux_nonconservative_ersing_etal)
surface_flux = (FluxHydrostaticReconstruction(FluxPlusDissipation(flux_ersing_etal,
                                                                  DissipationLocalLaxFriedrichs()),
                                              hydrostatic_reconstruction_ersing_etal),
                FluxHydrostaticReconstruction(flux_nonconservative_ersing_etal,
                                              hydrostatic_reconstruction_ersing_etal))

polydeg = 2
basis = LobattoLegendreBasis(polydeg)

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
# Get the TreeMesh and setup a non-periodic mesh

coordinates_min = -3.0
coordinates_max = 3.0
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 8,
                n_cells_max = 10000,
                periodicity = false)

# create the semidiscretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers

tspan = (0.0, 100.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false)

stepsize_callback = StepsizeCallback(cfl = 0.7)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 10,
                                     save_initial_solution = true,
                                     save_final_solution = true)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback,
                        stepsize_callback, save_solution)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(stage_limiter!), dt = 1.0,
            save_everystep = false, callback = callbacks, adaptive = false);

summary_callback() # print the timer summary
