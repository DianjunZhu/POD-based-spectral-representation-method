## POD-based-spectral-representation-method
A Matlab implementation of Proper Orthogonal Decomposition (POD) based spectral representation for simulating stochastic processes.

## Files Description

- **`POD_fft.m`** - Core implementation containing the main function for POD-based spectral representation simulation
- **`POD_1DStochasticProcess.m`** - Example 1: Simulation of non-stationary stochastic processes
- **`POD_1DGuassMotion.m`** - Example 2: Simulation of non-stationary seismic ground motions

## Quick Start

Run the examples directly

## Multi-Dimensional Extension

The implementation can be **easily extended to multi-dimensional cases**:
for 2D case:
S_matrix = [[S_xx(ω,t), S_xy(ω,t)],
            [S_yx(ω,t), S_yy(ω,t)]]
The same POD framework applies directly according to Ref:
Huang, G. 2015. Application of proper orthogonal decomposition in fast Fourier transform—Assisted multivariate nonstationary process simulation.

