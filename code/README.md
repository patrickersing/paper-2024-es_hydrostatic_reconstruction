# code

This folder contains code to reproduce the results in Section 4 of the article "Entropy-stable hydrostatic reconstruction schemes for shallow water systems".
It contains:
- The `Manifest.toml` and `Project.toml` of the Julia project.
- A modified version of the unregistered `TrixiShallowWater.jl` package.
- The file `noncons_bc.jl` with modifications to enable nonconservative boundary conditions on unstructured meshes in `Trixi.jl`.
- The sub-directories `Convergence_ML`, `DamBreakTriangular`, `ML_Dambreak`, `WB_Curvilinear`, `Lock_exchange`, `Parabolic_bowl_2d`, and `WB_ML` that contain setup files and visualization routines for the numerical experiments.
