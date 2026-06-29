# Intraseasonal-Planetary-Equatorial-Synoptic-Dynamics-IPESD-model

## Overview

This project provides MATLAB scripts for simulating the **Intraseasonal Planetary Equatorial Synoptic Dynamics (IPESD) model**, a multi-scale model designed to study the upscale effects of **Convectively Coupled Kelvin Waves (CCKWs)** on the Madden-Julian Oscillation (MJO). The method follows the multi-scale governing equations derived by **Majda & Biello (2004, 2005)**, coupling synoptic-scale perturbations with planetary-scale responses through **Equatorial Momentum Transport (EMT)** and **Equatorial Heat Transport (EHT)**.

Key features:

* Multi-scale decomposition of atmospheric variables: weather-scale perturbations and planetary-scale envelopes.
* Explicit representation of vertical heating profiles for deep convection and shallow clouds.
* Analytic computation of EMT and EHT as forcing terms for planetary-scale waves.
* Quasi-linear solution method for efficient simulation of large-scale MJO responses.
* Adjustable synoptic-scale heating profiles to explore different CCKW structures.

---

## Model Setup

1. **Synoptic-scale heating — governed by the synoptic-scale equatorial weak temperature gradient (SEWTG) equations:**  
   * Defines vertical and horizontal heating profiles for deep convection and shallow clouds.
   * Heating region is idealized in the Indian Ocean with a width of 5000 km.
   * Includes first and second baroclinic modes representing deep and shallow convection.

2. **Planetary-scale response — governed by the quasilinear equatorial long-wave (QLELWE) equations:**  
   * Driven by EMT and EHT computed from SEWTG solutions.
   * Solved using spectral Hermite modal decomposition for baroclinic modes.

3. **Coupling Mechanism:**  
   * Synoptic heating → EMT/EHT → forcing of planetary-scale equations → MJO response.

---

## Scripts Overview

| Script                     | Description                                                                                   |
| --------------------------- | --------------------------------------------------------------------------------------------- |
| `IPESD-syn.m`              | Defines synoptic-scale heating; computes velocity, temperature, and pressure perturbations.   |
| `IPESD-main.m`             | Solves planetary-scale QLELWE using Hermite spectral decomposition; generates response fields. |
| `data/`                     | Stores intermediate and output matrices (`l`, `v`, `r`) representing planetary-scale response. |

---

## Example Usage

```matlab
% Step 1: Define synoptic-scale heating
run('IPESD-syn.m');

% Step 2: Solve planetary-scale equations
run('IPESD-main.m');

% Step 3: Load results and visualize
load('data/l01.txt'); % Example output: zonal velocity
imagesc(l); colorbar; xlabel('Longitude'); ylabel('Time');
```

## Users Guide

Users can modify synoptic heating amplitudes and phases to simulate different CCKW structures.  
Output matrices represent planetary-scale responses suitable for further analysis.

---

## Key Parameters

| Parameter        | Description               | Value / Unit          |
|-----------------|---------------------------|---------------------|
| e               | Froude number             | 0.125               |
| Lx, Ly          | Warm pool dimensions       | 5000 km × 2000 km  |
| z_max           | Tropopause height         | 16 km               |
| u_scale         | Horizontal velocity scale | 6.25 m/s            |
| w_scale         | Vertical velocity scale   | 0.025 m/s           |
| momentum_damping| Momentum damping          | 0.18 day⁻¹          |
| heat_damping    | Thermal damping           | 0.1 day⁻¹           |

*(Complete parameters are consistent with Biello & Majda 2005.)*

---

## Figures

### Synoptic-scale Heating and Velocity Fields
<p align="center">
  <img src="figures/synoptic_heating_velocity.jpg" width="70%" />
</p>

### EMT and EHT
<p align="center">
  <img src="figures/emt_eht.jpg" width="80%" />
</p>

### Planetary-scale Response: Fast
<p align="center">
  <img src="figures/Planetary-Response-fast.png" width="90%" />
</p>

### Planetary-scale Response: SLow
<p align="center">
  <img src="figures/Planetary-Response-slow.png" width="90%" />
</p>

---

## References

* Biello, J. A., and A. J. Majda, 2005: A New Multiscale Model for the Madden–Julian Oscillation. Journal of the Atmospheric Sciences, 62, 1694–1721, https://doi.org/10.1175/JAS3455.1. *
* Majda, A. J., and J. A. Biello, 2004: A multiscale model for tropical intraseasonal oscillations. Proceedings of the National Academy of Sciences, 101, 4736–4741, https://doi.org/10.1073/pnas.0401034101. *

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

