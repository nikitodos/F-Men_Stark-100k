# STARK-100k — Solid Rocket Motor Design & Analysis

<p align="center">
  <img src="Logo_Project.png" width="600">
</p>

Politecnico di Milano · Space Propulsion · A.Y. 2025/26  
**F-Men: Days of Future Blast**

---

## Overview

This repository contains the MATLAB toolchain developed for the parametric design and analysis of **STARK-100k**, a 100 kN-class AP/HTPB solid rocket motor.

The workflow combines internal ballistics, grain sizing, nozzle flow and thermal analysis, active water cooling, feed-system sizing, and uncertainty analysis into a modular engineering workflow.

### Key features

- AP/HTPB propellant ballistic fitting using the Saint-Robert law
- Parametric BATES grain design
- Quasi-steady firing simulation
- Isentropic nozzle and area–Mach reconstruction
- Transient nozzle thermal analysis
- Bartz heat-transfer correlation with curvature correction
- TBC and MDN250 thermal modelling
- Active water-cooling and pressure-drop analysis
- Feed-line sizing and water-hammer analysis
- Monte Carlo uncertainty analysis
- 2D axisymmetric motor visualisation

## Nominal Design

| Parameter | Value |
|---|---:|
| Propellant | AP/HTPB 78/22 |
| Sea-level thrust | 100 kN |
| Chamber pressure | 70 bar |
| Throat diameter | 107 mm |
| Chamber radius | 186 mm |
| Grain length | 1187 mm |
| Burn time | 24.6 s |
| Delivered Isp | 245 s |
| MEOP, 99th percentile | 74.7 bar |
| Maximum MDN250 temperature | 729 K |
| Coolant mass flow | 3.37 kg/s |
| Total wet cooling-system mass | 188 kg |

## Repository Structure

```text
STARK-100k/
├── main.m
├── Firing.m
├── SolidUnitDesign.m
├── Nozzle_Transient_Heatmap.m
├── feed_line_design.m
├── refine_cooling_support.m
├── run_MC.m
├── uncertainty.m
├── render_engine_section_2D.m
├── [supporting MATLAB functions]
├── Stark-100k_FINAL_PLOT.png
└── FMen_DAloisio_et_al_2026.pdf
```

## Getting Started

### Requirements

- MATLAB R2025b or later
- No additional toolbox is required for the core workflow

### Run

```matlab
main
```

The main script performs the integrated design workflow and prints the principal results to the MATLAB command window.

## Academic Context

**Course:** Space Propulsion  
**Institution:** Politecnico di Milano  
**Academic Year:** 2025/2026

### Team

F-Men: Days of Future Blast

- Giovanni Nicola D'Aloisio
- Marco Gattone
- Nikolaas Valentin Pacheco
- Simone Rossi
- Tommaso Elia Tempesta

## Report

The complete technical report is included in the repository:

`FMen_DAloisio_et_al_2026.pdf`

## Disclaimer

This repository contains academic engineering work. The models, assumptions, correlations, and results are intended for educational and preliminary engineering-analysis purposes and are not flight-qualified propulsion software.

## License & Attribution

The repository and associated technical material remain subject to the authorship and intellectual-property conditions stated in the accompanying report.
