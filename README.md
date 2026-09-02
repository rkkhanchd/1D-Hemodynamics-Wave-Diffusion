# 1D Hemodynamics: Wave-Diffusion Model

MATLAB implementation of a one-dimensional wave-diffusion model for blood-flow propagation through a healthy artery and an arterial segment containing an aneurysm.

The simulations investigate the effect of arterial wall viscoelasticity on wave propagation, attenuation, reflection, transmission, and delay.

## Overview

The numerical model solves the one-dimensional flow-rate perturbation equation using a smooth-coefficient physiological formulation. The coefficients vary smoothly between the healthy arterial region and the aneurysm region.

Three cases are simulated using the same:

- total arterial length
- inlet pressure pulse
- boundary conditions
- numerical method
- spatial mesh
- time step

Cases 2 and 3 use identical aneurysm geometry and wall stiffness/thickness parameters. They differ only in wall viscosity, allowing the effect of wall viscoelasticity to be isolated.

### Simulation cases

| Case | Arterial configuration | Wall model |
|---|---|---|
| Case 1 | Fully healthy artery | Elastic |
| Case 2 | Healthy–aneurysm–healthy | Elastic aneurysm |
| Case 3 | Healthy–aneurysm–healthy | Viscoelastic aneurysm |

## Model

The dimensionless flow-rate perturbation is defined as

$$
q=\frac{\tilde{Q}}{Q_\star},
$$

with dimensionless axial coordinate

$$
x=\frac{z}{L},
$$

and dimensionless time

$$
\tau=\frac{t}{T}.
$$

The governing wave-diffusion equation implemented in the numerical model is

$$
q_{\tau\tau}
+G(x)q_\tau
-C(x)^2q_{xx}
-V(x)q_{\tau xx}
=0.
$$

The model coefficients are obtained from the local arterial wall and reference-area properties.

The numerical implementation uses smoothly varying coefficients across the aneurysm necks. No internal jump conditions are imposed in the physiological numerical simulations.

## Geometry and material parameters

### Healthy arterial region

The healthy arterial radius, wall thickness, and elastic modulus are defined in the MATLAB implementation.

The healthy elastic modulus is calculated using the Olufsen empirical relation implemented in the code.

### Aneurysm region

The aneurysm parameters used in the simulations are:

- Radius: 2.50 mm
- Wall thickness: 0.15 mm
- Elastic modulus: 1.00 MPa
- Wall viscosity:
  - Case 2: 0
  - Case 3: 2.00 × 10³ Pa·s

The aneurysm occupies

$$
z_{a1}=72.5\,\mathrm{mm}
$$

to

$$
z_{a2}=77.5\,\mathrm{mm},
$$

within a total arterial length of 150 mm.

The transition between healthy and aneurysmal properties is implemented using a smooth localization function.

## Numerical setup

The code was developed and tested using **MATLAB R2024b**.

The simulations use the following common numerical parameters:

| Parameter | Value |
|---|---:|
| Total arterial length, $L$ | 0.150 m |
| Characteristic time, $T$ | 0.050 s |
| Blood density, $\rho$ | 1056 kg/m³ |
| Dynamic viscosity, $\mu$ | $3.5 \times 10^{-3}$ Pa·s |
| Requested spatial intervals | 500 |
| CFL safety factor | 0.9 |
| Final physical simulation time | 2.0 s |
| Requested physical time step | $2.0 \times 10^{-5}$ s |

The time step is selected using the CFL constraint and is shared by all three simulations.

## Inlet function

The inlet pressure perturbation is defined using a triple-sech pulse with a smooth ($C^\infty$) activation function.

The raw inlet pressure pulse is defined as

$$
p_{\mathrm{raw}}(t) =
\frac{P_P}{\cosh\left(\frac{t-t_{\mathrm{delay}}-t_P}{w_P}\right)}
+
\frac{P_T}{\cosh\left(\frac{t-t_{\mathrm{delay}}-t_T}{w_T}\right)}
+
\frac{P_D}{\cosh\left(\frac{t-t_{\mathrm{delay}}-t_D}{w_D}\right)}.
$$

The applied inlet pressure is

```math
p_{\mathrm{in}}(t) = S_\infty(t)\,p_{\mathrm{raw}}(t)
```

where $S_\infty(t)$ is the smooth activation function.

The dimensionless inlet boundary condition is

```math
f_2(\tau) = \frac{C_h\,p_{\mathrm{in}}(T\tau)}{\alpha_h A_\star}
```

The corresponding inlet flow-rate perturbation is

```math
\tilde{Q}_{\mathrm{in}}(t) = \frac{c_h}{\alpha_h}\,p_{\mathrm{in}}(t)
```
The inlet-function check produces:

1. $p_{\mathrm{in}}(T\tau)$
2. $\tilde{Q}_{\mathrm{in}}(T\tau)$
3. $f_2(\tau)$
4. $S_\infty(t)$

The code also reports:

- $\max p_{\mathrm{in}}$
- $\max |\tilde{Q}_{\mathrm{in}}|$
- $\max |f_2|$
- dimensional and dimensionless times at which the maxima occur
- $f_2(0)$
- numerical estimates of $f_2'(0)$

## Repository structure

1D-Hemodynamics-Wave-Diffusion/
│
├── README.md
│
├── code/
│   └── healthy_vs_aneurysm_comparison.m
│
├── results/
│   ├── inlet_function_check.png
│   ├── case1_full_healthy_probe_waveforms.png
│   ├── case1_full_healthy_space_time_maps.png
│   ├── case2_elastic_aneurysm_probe_waveforms.png
│   ├── case2_elastic_aneurysm_space_time_maps.png
│   ├── case3_viscoelastic_aneurysm_probe_waveforms.png
│   └── case3_viscoelastic_aneurysm_space_time_maps.png
│
└── .gitignore
