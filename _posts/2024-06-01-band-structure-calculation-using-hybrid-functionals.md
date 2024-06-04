---
title: 'Band-structure calculation using hybrid functionals'
date: 2024-06-01
tags:
  - VASP
  - band structure
  - hybrid functional
---

# [Band-structure calculation using hybrid functionals](https://www.vasp.at/wiki/index.php/Band-structure_calculation_using_hybrid_functionals)

> **a significant difference between hybrid band-structure calculations and DFT band-structure calculations**:
> 
> The electronic charge density suffices for density functionals to define the Hamiltonian, and no regular **k** mesh is required during DFT band-structure calculations. However, if no regular **k** mesh is provided, the electronic charge density must be fixed during the DFT band-structure calculation by setting [ICHARG](https://www.vasp.at/wiki/index.php/ICHARG "ICHARG")=11 in the [INCAR](https://www.vasp.at/wiki/index.php/INCAR "INCAR") file.
> 
> **Warning:** The electronic charge density must not be fixed for any hybrid calculation, i.e., never set [ICHARG](https://www.vasp.at/wiki/index.php/ICHARG "ICHARG")=11!

根据以上信息，采用杂化泛函不能使用ICHARG=11，即只能做自洽计算。