---
title: 'Perovskite Framework Transformations'
date: 2025-3-6
tags:

  - Perovskite
---

# Perovskite Framework Transformations

## Basics

[***Aristotype***](https://dictionary.iucr.org/Aristotype) - Cubic $Pm\bar3m$ (221)

> An **aristotype** is a high-symmetry structure type that can be viewed as an idealized version of a lower-symmetry structure.
>
> The lower-symmetry structure is called **hettotype**.
>
> Aristotypes are also known as *basic structures* and hettotypes as *derivative structures*.

- Corner-connected $BX_6$ octahedral framework

- Face-sharing $AX_{12}$ cuboctahedra

**$ABX_3$ Perovskite Framework**

| Cubic Close Packed (CCP) | Hexagonal Close Packed (HCP) |
| ------------------------ | ---------------------------- |
| ...C(ABC)A...            | ...B(AB)A...                 |
| $Pm\bar3m$ 221           | $r\bar3c$ 167                |
| **Aristotype**           | **Hettotype**                |

Topological transformation mechanism: **Anti-parallel $<111>_\text{Cubic}$ Rotation**

## Quantitative Analysis of $BX_6$ Tilting

Glazer Notation:

Untilted - $a^0/b^0/c^0$

Tilted in the same direction - $a^+/b^+/c^+$

Tilted in opposite directions - $a^-/b^-/c^-$

***Ref:***

[The classification of tilted octahedra in perovskites - Glazer - 1972 - Acta Crystallographica Section B - Wiley Online Library](https://onlinelibrary-wiley-com.remotexs.ntu.edu.sg/doi/abs/10.1107/S0567740872007976)

[(IUCr) Some structures topologically related to cubic perovskite (E21), ReO3 (D09) and Cu3Au (L12)](https://journals.iucr.org/paper?S0567740877012114)

## Cubic $Pm\bar3m$ (221) Aristotype Framework

**Putative Composition**: $ABX_3$

**Unit Cell Parameter**: **a** = 2*d*, where *d* is the octahedral B-X bond length

**Unit Cell Volume**: ***V***(0) = $(2d)^3$

**Unit Cell Contents**:

| Species | Wyckoff symbol | coordinates                     |
| ------- | -------------- | ------------------------------- |
| A       | *1b*           | $\frac12$, $\frac12$, $\frac12$ |
| B       | *1a*           | 0, 0, 0                         |
| X       | *3d*           | $\frac12$, 0, 0                 |

**Ideal Kernel:** $AX_{12}$ **cuboctahedra** - 2 $\times$ hemi-cuboctahedron ($AX_9$)

> Rotating the upper part of a **twin-cuboctahedron** with a mirror plane generates a **hemi-cuboctahedron**, in which case the mirror plane is destroyed.

## O’Keeffe & Hyde Tilt System Approach

- Rotation axes are specified for the eight octahedra with centers at the corners of the perovskite unit cell, i.e. for a 2 $\times$ 2 $\times$ 2 'supercell' in which eight octahedra surrounding an $AX_{12}$ coordination sphere.
- Rotation is clockwise when viewed in the direction of the axis, e.g. [111] and [$\bar1\bar1\bar1$] represent rotations about the same axis, but in opposite senses.
- The size of the octahedra is maintained constant, the B-X distance being *d* and the octahedron edge length $\sqrt2d$

Four rotation systems are considered:

1. Trigonal $R\bar3c$ (167)
2. Cubic $Im\bar3$ (204)
3. Tetragonal $I4/mmm$ (123)
4. Orthorhombic $Pnma$ (62)

Each of these different tilt systems gives rise to different ***$AX_{12}$ coordination spheres***, which are considered to be the kernels of the perovskites.

## Perovskite Hettotypes

### 1. Trigonal $R\bar3c$ (167)

**Putative Composition**: $ABX_3$

**Rotation System**: $<111>$

**Unit Cell Parameter (hexagonal metric)**:

**a** = $\sqrt8d\cos\phi$

**c** = $\sqrt{48}d$

where *d* is the octahedral edge length

**Unit Cell Volume**: ***V***($\phi$) = ***V***(0)$\cos^2\phi$

**Unit Cell Contents**:

| Species | Wyckoff symbol | coordinates                                                  |
| ------- | -------------- | ------------------------------------------------------------ |
| A       | *6a*           | 0, 0, $\frac14$                                              |
| B       | *6b*           | 0, 0, 0                                                      |
| X       | *18e*          | $x$, 0, $\frac14$ with $x = \frac{(\sqrt3-\tan\phi)}{\sqrt{12}}$ |

**Ideal Kernel:** $AX_{12}$  **[triangular bifrustum](https://en.wikipedia.org/wiki/Triangular_bifrustum)** - 2 $\times$ triangular frustum

**Topological Transformation**: $Pm\bar3m$ (221) $\rightarrow$ $R\bar3c$ (167)

Anti-parallel rotation - $a^-$

At $\phi=30\degree$ conversion from CCP to HCP complete with

| Species | Wyckoff symbol | coordinates             |
| ------- | -------------- | ----------------------- |
| X       | *18e*          | $\frac13$, 0, $\frac14$ |

### 2. Cubic $Im\bar3$ (204)

**Putative Composition**: $A^\prime A^{\prime\prime}_3B_4X_{12}$

**Rotation System**: 4 $\times$ $<111>$ directions

**Unit Cell Parameter**:

**a** = $\frac{8\cos\phi+4}3d$, where *d* is the octahedral edge length

**Unit Cell Volume**: ***V***($\phi$) = $\frac{(2\cos\phi+1)^3}{27}$***V***(0)

**Unit Cell Contents**:

| Species            | Wyckoff symbol | coordinates                                                  |
| ------------------ | -------------- | ------------------------------------------------------------ |
| $A^\prime$         | *2a*           | 0, 0, 0                                                      |
| $A^{\prime\prime}$ | *6b*           | 0, $\frac12$, $\frac12$                                      |
| B                  | *8c*           | $\frac14$, $\frac14$, $\frac14$                              |
| X                  | *24g*          | 0, $y$, $z$ with $y = \frac{3\cos\phi+\sqrt3\sin\phi}{8\cos\phi+4}$ $z = \frac{3\cos\phi-\sqrt3\sin\phi}{8\cos\phi+4}$ |

**Ideal Kernel:** $AX_{12}$  **[icosahedron](https://en.wikipedia.org/wiki/Icosahedron)**

**Topological Transformation**: $Pm\bar3m$ (221) $\rightarrow$ $Im\bar3$ (204)

Parallel rotation - $a^+$

At $\phi=22.24\degree$ there is perfect $AX_{12}$ icosahedron with

| Species | Wyckoff symbol | coordinates     |
| ------- | -------------- | --------------- |
| X       | *24g*          | 0, 0.301, 0.186 |

### 3. Tetragonal $I4/mmm$ (123)

**Putative Composition**: $A^\prime A^{\prime\prime} A^{\prime\prime\prime}_2B_4X_{12}$

**Rotation System**: $<110>$

**Unit Cell Parameter**:

**a** = $2d(1+\cos\phi)$

**c** = $4d\cos\phi$

where *d* is the octahedral edge length

**Unit Cell Volume**: ***V***($\phi$) = $\frac{(1+\cos^2\phi)}4\cos\phi$***V***(0)

**Unit Cell Contents**:

| Species                  | Wyckoff symbol | coordinates                                                  |
| ------------------------ | -------------- | ------------------------------------------------------------ |
| $A^\prime$               | *2a*           | 0, 0, 0                                                      |
| $A^{\prime\prime}$       | *2b*           | 0, 0, $\frac12$                                              |
| $A^{\prime\prime\prime}$ | *4c*           | 0, $\frac12$, 0                                              |
| B                        | *8f*           | $\frac14$, $\frac14$, $\frac14$                              |
| $X^\prime$               | *8h*           | $x$, $x$, 0 with $x = \frac{1+\cos\phi+\sqrt2\sin\phi}{4+4\cos\phi}$ |
| $X^{\prime\prime}$       | *16n*          | 0, $y$, $z$ with $y = \frac1{2+2\cos\phi}$ $z = \frac{\sqrt2-\tan\phi}{\frac4{\sqrt2}}$ |

**Ideal Kernel:**

- $AX_{12}$ face-sharing [square antiprism](https://en.wikipedia.org/wiki/Square_antiprism): 2 $\times$ square antiprism ($A^{\prime\prime}X_8$)
- $AX_{12}$ tetracapped cube: cube ($A^\prime X_8$) + 4 $\times$ tetrahedron

**Topological Transformation**: $Pm\bar3m$ (221) $\rightarrow$ $I4/mmm$ (123)

Parallel rotation - $a^+/b^+$

At $\phi=19.47\degree$ there is near perfect $AX_{12}$  face-sharing square antiprisms and tetracapped cubes with

| Species            | Wyckoff symbol | coordinates       |
| ------------------ | -------------- | ----------------- |
| $X^\prime$         | *8h*           | 0.3107, 0.3107, 0 |
| $X^{\prime\prime}$ | *16n*          | 0, 0.2574, 0.1875 |

### 4. Orthorhombic $Pnma$ (62)

**Putative Composition**: $A_4B_4X^\prime_4X^{\prime\prime}_8$

**Rotation System**: $<0\bar11>$

**Unit Cell Parameter**:

**a** = $d[\frac{8(2+\cos^2\phi)}3]^{\frac12}$

**b** = $d[\frac{48}{1+\sec^2\phi}]^{\frac12}$

**c** = $d\sqrt8\cos\phi$

where *d* is the octahedral edge length

**Unit Cell Volume**: ***V***($\phi$) = $\cos^2\phi$***V***(0)

**Unit Cell Contents**:

| Species            | Wyckoff symbol | coordinates                                                  |
| ------------------ | -------------- | ------------------------------------------------------------ |
| A                  | *4c*           | $x$, $\frac14$, $z$                                          |
| B                  | *4b*           | 0, 0, $\frac12$                                              |
| $X^\prime$         | *4c*           | $x$, $x$, 0 with $x = \frac{(\cos^2\phi-1)}{2\cos^2\phi+4}$ $z = \frac{\sqrt3+\tan\phi}{\sqrt{12}}$ |
| $X^{\prime\prime}$ | *8d*           | $x$, $y$, $z$ with $x = \frac{2-\sqrt3\sin\phi\cos\phi+\cos^2\phi}{8+4\cos^2\phi}$ $y = -\frac{\tan\phi}{\sqrt{48}}$ $z = \frac{3\sqrt3+\tan\phi}{\sqrt{48}}$ |

**Ideal Kernel:** $AX_{12}$ augmented tetracapped trigonal prism (approx)

$AX_6$ trigonal prism $\rightarrow$ $AX_9$ tricapped trigonal prism $\rightarrow$ $AX_{10}$ tetracapped trigonal prism $\rightarrow$ $AX_{12}$ augmented tetracapped trigonal prism

**Topological Transformation**: $Pm\bar3m$ (221) $\rightarrow$ $Pnma$ (62)

Anti-parallel rotation - $a^-/b^+$

For $\phi=30\degree$ the $BX_6$ octahedra are perfectly regular and can be described as twinned HCP. Trigonal prisms ($AX_6$) are created between the twin planes.

> **Anti-perovskite**
>
> In cementite $Fe_3C$ only the $CFe_6$ trigonal prisms are observed at the twin planes and the octahedral sites are unoccupied.

## Tolerance Factor

A tolerance factor ($t_p$) for perovskite defines the size constraints for A cations to fit exactly inside the cavities in the  $BX_6$ framework, i.e. $AX_{12}$ unit.
$$
(r_A+r_X) = \sqrt2(r_B+r_X)\\
t_p = \frac{r_A+r_X}{\sqrt2(r_B+r_X)}
$$
For a perfect fit $t_p=1.0$ but perovskites exist from ~0.8 - 1.0

Complications include:

- Coordination number and effective bond valence
- Size as function of temperature and pressure
- Molecular dynamics of organic component
- Hydrogen bonding

***Ref:***

[An extended Tolerance Factor approach for organic–inorganic perovskites - Chemical Science (RSC Publishing)](https://pubs.rsc.org/en/content/articlelanding/2015/sc/c5sc00961h)

## Summary

- The tilting of the octahedra leads to a densification of the structure.
- The patterns of $BX_6$ octahedral tilting in 3D $ABX_3$ perovskites will follow regular hierarchical sequences.
- While in some cases the $BX_6$ octahedra are regular (remain as rigid bodies), there is sometimes distortion required or a small secondary tilt used to maintain the corner-connected topology.
- The same rules apply to the slabs (even single octahedral sheets) of 'intercalated' perovskite that are 2D structures.

## Key Text Book

[Perovskites Modern and Ancient. By Roger H. Mitchell. Thunder Bay, Ontario: Almaz Press , 2002. Price USD 70.00. ISBN 0-9689411-0-9](https://onlinelibrary.wiley.com/iucr/doi/10.1107/S0108768102020220)

Perovskites and High Tc Superconductors. Francis S. Galasso. ISBN 978288124391

[Crystal Structures: A Working Approach. Helen D. Megaw: 9780721662602](https://www.amazon.com/Crystal-Structures-Working-Approach-chemistry/dp/0721662609)
