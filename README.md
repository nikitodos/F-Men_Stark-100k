# STARK-100k – Solid Rocket Motor Design & Analysis Suite

[![MATLAB](https://img.shields.io/badge/MATLAB-R2025b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

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

All analyses are integrated into a single main script (`main.m`) that produces **12 figures** and a detailed console summary.

---

## Repository structure
