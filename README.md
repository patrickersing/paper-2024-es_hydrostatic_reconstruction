# Entropy-stable hydrostatic reconstruction schemes for shallow water systems

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

This reproducibility repository contains the information and code to reproduce the results of the article 

**Placeholder: Paper**

If you find these results useful, please cite the article mentioned above. If you use the implementations provided here, please also cite this repository as

**Placeholder: Zenodo**

## Abstract
In this work, we develop a new hydrostatic reconstruction procedure  to construct well-balanced schemes for one and multilayer shallow water flows, including wetting and drying. Initially, we derive the method for a path-conservative finite volume scheme and combine it with entropy-conservative fluxes and suitable numerical dissipation to preserve an entropy inequality in the semi-discrete case. We then combine the novel hydrostatic reconstruction with a collocated nodal split-form discontinuous Galerkin spectral element method, extending the method to high-order and curvilinear meshes. The high-order method incorporates an additional positivity-limiter and is blended with a compatible subcell finite volume method to maintain well-balancedness at wet/dry fronts. We prove entropy-stability, well-balancedness and positivity-preservation for both methods.Numerical results for the high-order method validate the theoretical findings and demonstrate the robustness of the scheme.

## Installation
1. Install Julia v1.10
2. Download the repository
```
git clone https://github.com/patrickersing/paper-2024-es_hydrostatic_reconstruction_dev.git
```
1. Set the working directory to `/code` and instantiate the Julia environment using
```bash
cd paper-2024-es_hydrostatic_reconstruction/code
julia --project=. -e 'using Pkg; Pkg.develop(PackageSpec(path="TrixiShallowWater.jl")); Pkg.instantiate()'
```

## Usage
To reproduce the results in the paper start a new `Julia` session with the `--project=.` flag. It is possible (but optional) to use multiple threads by appending the `--threads=N` flag, to run on `N` threads. 
```bash
julia --project=. --threads=N
```

Then execute the following commands in the Julia-REPL to create the respective results:
- **Section 4.1 - Figure 4**
    The following command will create the figure and save it as a `.pdf` within the `/Convergence_ML` folder 
    ```julia
    cd("Convergence_ML/")
    include("visualization_spectral_convergence.jl")
    cd("../")
    ```

- **Section 4.2 - Figure 5 + 6**
    The following command will create Figure 5 + 6 and save them as a `.pdf` within the `WB_ML/` folder
    ```julia
    cd("WB_ML/")
    include("visualization_setup.jl")
    include("visualization_time_series.jl")
    cd("../")
    ```

- **Section 4.3 - Figure 7 + Table 1**
    The following command will create Figure 7 and save it as a `.pdf` within the `WB_Curvilinear/` folder. Results from Table 1 are displayed in the REPL. 
    ```julia
    cd("WB_Curvilinear/")
    include("visualization_bottom_topography.jl")
    cd("../")
    ```

- **Section 4.4 - Figure 8 + 9**
    The following command will create Figure 5 + 6 and save them as a `.pdf` within the `WB_ML/` folder
    ```julia
    cd("DamBreakTriangularBottom/")
    include("visualization_setup.jl")
    include("visualization_gauge_data.jl")
    cd("../")
    ```
  
- **Section 4.5 - Figure 10**
    The following commands creates output files for Figure 8 and converts them to `.vtu` format readable in Paraview
    ```julia
    cd("ML_Dambreak/")
    include("elixir_shallowwater_multilayer_dam_break_dry.jl")
    # convert .h5 output files to .vtu
    using Trixi2Vtk
    trixi2vtk("out/solution*", output_directory="out", reinterpolate=false)
    cd("../")
    ```
    To create Figure 8 then start `Paraview 5.10.0-RC1` and load the statefile `ML_Dambreak/visu.pvsm`

## Authors

- [Patrick Ersing](https://liu.se/en/employee/pater53) (Linköping University, Sweden)
- [Andrew R. Winters](https://liu.se/en/employee/andwi94) (Linköping University, Sweden)

## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
