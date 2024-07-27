---
title: 'VASP Error & Solution'
date: 2024-03-10
tags:
  - VASP
---

<div style="color:black; background-color:#CCCCFF; border: 1px solid #FFE0C3; border-radius: 10px; margin-bottom:0rem">
    <p style="margin:1rem; padding-left: 1rem; line-height: 2.5;">
    <b>Task Type:<\b> Elastic constant calculation
    <b>Error log:<\b> ` internal error in GENERATE_KPOINTS_TRANS: number of G-vector changed in star
        X        Y
`
    <b>Solution:<\b> Copy `IBZKPT` to `KPOINTS` and recalculate
    </p>
</div>