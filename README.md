# Offshore Wind Foundation – Demo Code

## Overview
This repository contains a **demonstration** of numerical simulation for multidirectional loading in offshore wind turbine foundations (monopod suction bucket).  
The workflow integrates **CPT-based soil resistance** and **soil–structure interaction models** to evaluate foundation performance under lateral and rotational loads.  

⚠️ **Note:** The full validated framework and parametric study results are currently under journal review and will be released upon acceptance.
This repo is intended only to showcase the modeling approach and code structure.  

- **apdefiner.m** – Calculates the coefficient *α* to account for the effect of lateral soil resistance on interface normal stress (nested in `tzcurve.m`).  
- **BETA13.m** – Calculates coefficients *β₁* and *β₃* for the p–y curve (nested in `pycurve.m`).  
- **BETA24.m** – Calculates coefficients *β₂* and *β₄* for the p–y curve (nested in `pycurve.m`).  
- **coeffinder.m** – Evaluates p–y spring forces.  
- **coeffinder_deadload.m** – Computes p–y forces affecting t–z forces during settlement (before rotation).  
- **fi_finder.m** – Determines the peak friction angle (as a function of soil density ratio).  
- **points.m** – Defines spring node locations along the bucket skirt, skirt tip, and lid.  
- **pycurve.m** – Applies the p–y curve theory (calls `BETA13.m` and `BETA24.m`).  
- **qzcurve2.m** – Calculates q–z spring forces at the skirt tip.  
- **qzcurve4.m** – Calculates q–z spring forces beneath the lid.  
- **Ranking.m** – Computes Rankine pressure (nested in `pycurve.m`).  
- **tzcurve.m** – Calculates t–z spring forces on the skirt.  
- **tzcurve_deadload.m** – Computes t–z forces during settlement (before rotation).  
- **transformation.m** – Computes positions of points in 3D Cartesian frame relative to the body frame (after rotation), based on the center of rotation *Q*.  
- **xtransformation.m** – Computes positions of points in 3D Cartesian frame relative to the fixed frame (after rotation), based on the center of rotation *Q*.  
- **Zvalues.m** – Calculates the depth of each spring embedded in the soil.  

---

## Usage   
1. Clone this repository or download the `.zip`.  
2. Open MATLAB and run `MAINCODE.m`.  
3. A sample load–displacement (H–θ / M–θ) curve will be generated using representative soil parameters, adapted from [Wang et al. (2018)](https://doi.org/10.1016/j.oceaneng.2017.12.006). 

---

## Disclaimer
This is a **demo version** of my research code.  
- No real parametric study or journal data is included.  
- All soil parameters are placeholders for demonstration only.  
- The complete validated framework will be published after journal acceptance.  

---

## Author
**Sina Hadadi**  
Ph.D. Candidate – Wind Energy (Offshore Foundations)  
Kunsan National University, South Korea  
Email: sinahdme@gmail.com  
GitHub: [github.com/sinahdme](https://github.com/sinahdme)  
