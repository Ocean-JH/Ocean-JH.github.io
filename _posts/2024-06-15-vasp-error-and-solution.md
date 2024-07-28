---
title: 'VASP Error & Solution'
date: 2024-06-15
tags:
  - VASP
---

    Task Type: Elastic constant calculation
    Error log:

> ` internal error in GENERATE_KPOINTS_TRANS: number of G-vector changed in star
>         X        Y`
    Solution: Copy `IBZKPT` to `KPOINTS` and recalculate

------