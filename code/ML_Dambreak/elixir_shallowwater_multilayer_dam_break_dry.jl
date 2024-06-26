
using OrdinaryDiffEq
using Trixi
using TrixiShallowWater

###############################################################################
# Semidiscretization of the multilayer shallow water equations for a dam break test over a dry domain
# with a discontinuous bottom topography function
equations = ShallowWaterMultiLayerEquations2D(gravity_constant = 9.81,
                                              rhos = (0.9, 0.95, 1.0),
                                              threshold_desingularization = 1e-6)

# This test case uses a special work around to setup a truly discontinuous bottom topography 
# function and initial condition. First, a dummy initial_condition_dam_break is introduced to create
# the semidiscretization. Then the initial condition is reset with the true discontinuous values 
# from initial_condition_discontinuous_dam_break.

function initial_condition_dam_break(x, t, equations::ShallowWaterMultiLayerEquations2D)
    # Bottom topography
    b = 1.4 * exp(-10.0 * ((x[1])^2 + (x[2])^2))

    if x[1] < 0.0
        H = [1.0, 0.8, 0.6]
    else
        H = [b, b, b]
    end

    v1 = zero(H)
    v2 = zero(H)

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

    return prim2cons(SVector(H..., v1..., v2..., b),
                     equations)
end

initial_condition = initial_condition_dam_break

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
# Get the unstructured quad mesh from a file (downloads the file if not available locally)

mesh_file = joinpath(@__DIR__, "WarpedMesh.mesh")

mesh = UnstructuredMesh2D(mesh_file, periodicity = false)

# Boundary conditions
boundary_condition = Dict(:Top => boundary_condition_slip_wall,
                          :Left => boundary_condition_slip_wall,
                          :Right => boundary_condition_slip_wall,
                          :Bottom => boundary_condition_slip_wall)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition,
                                    solver, boundary_conditions = boundary_condition)

###############################################################################
# ODE solver

tspan = (0.0, 2.0)
ode = semidiscretize(semi, tspan)

###############################################################################
# Workaround to set a discontinuous bottom topography and initial condition.

# alternative version of the initial condition used to setup a truly discontinuous
# test case and initial condition.
# In contrast to the usual signature of initial conditions, this one get passed the
# `element_id` explicitly. In particular, this initial conditions works as intended
# only for the specific mesh loaded above!
function initial_condition_discontinuous_dam_break(x, t, element_id,
                                                   equations::ShallowWaterMultiLayerEquations2D)
    # Bottom topography
    b = 1.4 * exp(-10.0 * ((x[1])^2 + (x[2])^2))

    if x[1] < 0.0
        H = [1.0, 0.8, 0.6]
    else
        H = [b, b, b]
    end

    N = 20  # mesh elements per dimension
    IDs_left = collect((N / 2):N:(N / 2 + (N - 1) * N))   # left element ids
    IDs_right = IDs_left .+ 1   # right element ids
    if element_id in IDs_left
        H = [1.0, 0.8, 0.6]
    elseif element_id in IDs_right
        H = [b, b, b]
    end

    # It is mandatory to shift the water level at dry areas to make sure the water height h
    # stays positive. The system would not be stable for h set to a hard 0 due to division by h in
    # the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    # with a default value of 5*eps() ≈ 1e-15 in double precision, is set in the constructor above
    # for the ShallowWaterMultiLayerEquations2D and added to the initial condition if h = 0.
    # This default value can be changed within the constructor call depending on the simulation setup.
    for i in reverse(eachlayer(equations))
        if i == nlayers(equations)
            H[i] = max(H[i], b + equations.threshold_limiter)
        else
            H[i] = max(H[i], H[i + 1] + equations.threshold_limiter)
        end
    end

    v1 = zero(H)
    v2 = zero(H)

    return prim2cons(SVector(H..., v1..., v2..., b),
                     equations)
end

# point to the data we want to augment
u = Trixi.wrap_array(ode.u0, semi)
# reset the initial condition
for element in eachelement(semi.solver, semi.cache)
    for j in eachnode(semi.solver), i in eachnode(semi.solver)
        x_node = Trixi.get_node_coords(semi.cache.elements.node_coordinates, equations,
                                       semi.solver, i, j, element)
        u_node = initial_condition_discontinuous_dam_break(x_node, first(tspan), element,
                                                           equations)
        Trixi.set_node_vars!(u, u_node, equations, semi.solver, i, j, element)
    end
end

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 500
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = true,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (energy_total,
                                                                 energy_kinetic,
                                                                 energy_internal))

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 0.5,
                                     save_initial_solution = true,
                                     save_final_solution = true, output_directory = "out")

stepsize_callback = StepsizeCallback(cfl = 0.9)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback, save_solution,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (Trixi.waterheight,))

###############################################################################
# run the simulation

sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()..., callback = callbacks, adaptive = false, dt = 1.0);
summary_callback() # print the timer summary
