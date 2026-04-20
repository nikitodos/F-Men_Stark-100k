# STARK-100k – Solid Rocket Motor Design & Analysis Suite

[![MATLAB](https://img.shields.io/badge/MATLAB-R2025b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)

<p align="center">
<img src="Logo_Project.png" alt="Project Logo" width="250"/>
</p>

**Politecnico di Milano – Space Propulsion A.Y. 2025/26**  
*Team F-Men: Days of Future Blast*  
D’Aloisio · Gattone · Pacheco · Rossi · Tempesta

---

## Overview

This repository contains a complete MATLAB toolchain for the **parametric design, transient simulation, thermal analysis, and feed‑system sizing** of a 100 kN thrust‑class solid rocket motor. The code was developed as a course project and is structured to be modular, reproducible, and easy to extend.

### Key features

- **Saint‑Robert (Vieille) ballistic fitting** from experimental data (AP/HTPB)
- **CEA‑based gas properties** for six AP/HTPB compositions (75/25 … 80/20)
- **Optimal grain geometry solver** (cylindrical BATES grain, neutral burning)
- **Quasi‑steady firing simulation** (mass balance, pressure/time history, Isp, c*, thrust)
- **Nozzle isentropic flow reconstruction** (area‑Mach solver)
- **Transient thermal analysis of the nozzle** (NTTA):
  - Gas‑side Bartz HTC with curvature correction
  - 1D lumped‑capacitance model (TBC + MDN250 liner)
  - External natural convection (Churchill‑Chu) and radiation
- **Cooling jacket energy balance** (iterative mass flow, pressure drop, ONB criterion)
- **Refined feed system design** (ASME B31.3 pipe sizing, water hammer analysis)
- **Uncertainty quantification** (Monte Carlo with log‑normal a and normal n)
- **2D axisymmetric renderings** of the motor section

All analyses are integrated into a single main script (`main.m`) that produces several figures and a detailed console summary.

---

## Repository structure

```
STARK-100k-SRM/
├── main.m                         # top‑level script (run this)
├── Firing.m                       # quasi‑steady 0‑D transient simulation
├── Nozzle_Transient_Heatmap.m     # transient thermal analysis (nozzle)
├── SolidUnitDesign.m              # parametric design loop
├── solve_geometry_fzero2.m        # BATES grain geometry solver
├── mach_from_area.m               # area‑Mach bisection solver
├── bartz_curv.m                   # Bartz HTC with curvature correction
├── sigma_bartz.m                  # Bartz sigma correction factor
├── refine_cooling_support.m       # refined pressure cascade, tank sizing, mass roll‑up
├── feed_line_design.m             # pipe diameter, wall thickness, water hammer (ASME B31.3)
├── darcy_friction.m               # Darcy‑Weisbach friction (laminar/Colebrook)
├── thermal_chain.m                # (not used in final version – legacy)
├── render_engine_section_2D.m     # 2D axisymmetric cross‑section plot
├── uncertainty.m                  # linear‑regression uncertainty for a and n
├── run_MC.m                       # helper for Monte Carlo calibration
├── ternary_str.m                  # inline conditional strings
├── P_sat_water.m                  # Antoine equation for water vapour pressure
├── density_water.m                # water density (polynomial)
├── viscosity_water.m              # water dynamic viscosity
├── conductivity_water.m           # water thermal conductivity
├── cp_water.m                     # water specific heat
├── README.md                      # this file
└── FMen_DAloisio_et_al_2026.pdf   # full project report (detailed methodology & results)
```

---

## Getting Started

### Prerequisites

- MATLAB R2025b or later (older versions may work, but not tested)
- No additional toolboxes required – all code uses base MATLAB functions.
- (Optional) The Utopia font for report‑style figures – otherwise MATLAB defaults are used.

### Running the full analysis

1. Clone the repository.
2. Open MATLAB in the repository folder.
3. Run `main.m`.

The script will:
- Load the experimental burning‑rate data and fit the Vieille-Saint Robert's (VSR) law.
- Evaluate six AP/HTPB compositions and select the best performer (by a weighted score).
- Generate the grain geometry, run the firing simulation, and perform the nozzle thermal analysis.
- Size the cooling jacket, the water tank, the pressurant N₂ tank, and the feed pipe (ASME B31.3).
- Produce all figures and print a detailed summary in the command window.

Typical execution time is 10 minutes (with Monte Carlo sampling number N = 2000).

### Customisation

- To change the thrust or total impulse, edit `Thrust_real` and `I_tot` in `main.m`.
- To modify the propellant composition range, edit the vectors `mf_AP` and `mf_HTPB`.
- To use a different TBC or case material, modify the `Mat` struct in `main.m` and the corresponding properties in `Nozzle_Transient_Heatmap.m`.

---

## Results summary (nominal design)

| Parameter               | Value                  |
|------------------------|------------------------|
| Propellant             | AP/HTPB 78/22          |
| Thrust (sea level)     | 100 kN                 |
| Chamber pressure       | 70 bar                 |
| Throat diameter        | 107 mm                 |
| Chamber radius         | 186 mm                 |
| Grain length           | 1187 mm                |
| Burn time              | 24.6 s                 |
| Isp (delivered)        | 245 s                  |
| MEOP (99% percentile)  | 74.7 bar               |
| Max metal temperature (MDN250) | 729 K (456 °C) – below limit |
| Coolant water flow     | 3.37 kg/s              |
| Total wet system mass  | 188 kg (water + tanks) |

*Full results are printed in the console when `main.m` is executed.*

---

## References (key sources)

6. **Sutton & Biblarz (2016)** – onset of nucleate boiling (ONB) criterion  
2. **Bartz (1957)** – convective heat transfer in rocket nozzles  
3. **Churchill & Chu (1975)** – natural convection from horizontal cylinders  
4. **ASME B31.3** – process piping (wall thickness, allowable stresses)  
5. **Moss & Basic (2013)** – pressure vessel design manual  

---

## License

This repository and all associated content, including but not limited to text, figures, tables, diagrams, algorithms, and technical data, are the intellectual property of the Design Team. Each section of this work is the intellectual property of its respective lead author(s), with contributions from team members as specified in the Authorship Declaration included in the final report. For permission requests or inquiries regarding commercial use, please contact the Team members.

---

## Authors

**F-Men: Days of Future Blast** – Space Propulsion course, Politecnico di Milano  
- Giovanni Nicola D’Aloisio  
- Marco Gattone  
- Nikolaas Valentin Pacheco  
- Simone Rossi  
- Tommaso Elia Tempesta  

For any questions or suggestions, please open an issue or contact the Team.

---

*Last updated: April 2026*
