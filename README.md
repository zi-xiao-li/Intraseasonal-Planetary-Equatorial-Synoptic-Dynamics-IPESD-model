# Intraseasonal-Planetary-Equatorial-Synoptic-Dynamics-IPESD-model

## Overview

This project provides MATLAB scripts for simulating the **Intraseasonal Planetary Equatorial Synoptic Dynamics (IPESD) model**, a multi-scale model designed to study the upscale effects of **Convectively Coupled Kelvin Waves (CCKWs)** on the Madden-Julian Oscillation (MJO). The method follows the multi-scale equations derived by **Majda & Biello (2004, 2005)**, coupling synoptic-scale perturbations with planetary-scale responses through **Equatorial Momentum Transport (EMT)** and **Equatorial Heat Transport (EHT)**.

Key features:

* Multi-scale decomposition of atmospheric variables: weather-scale perturbations and planetary-scale envelopes.
* Explicit representation of vertical heating profiles for deep convection and shallow clouds.
* Analytic computation of EMT and EHT as forcing terms for planetary-scale waves.
* Quasi-linear solution method for efficient simulation of large-scale MJO responses.
* Adjustable synoptic-scale heating profiles to explore different CCKW structures.

---

## Model Setup

1. **Synoptic-scale heating (SEWTG):**  
   * Defines vertical and horizontal heating profiles for deep convection and shallow clouds.
   * Heating region is idealized in the Indian Ocean with a width of 5000 km.
   * Includes first and second baroclinic modes representing deep and shallow convection.

2. **Planetary-scale response (QLELWE):**  
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
| `utils/`                    | Helper functions for Hermite modes and numerical operations.                                  |

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
