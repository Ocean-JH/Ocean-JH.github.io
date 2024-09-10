---
title: 'Ab-initio Molecular Dynamics （VASP）'
date: 2024-9-9
tags:
  - VASP
  - Molecular Dynamics
---

# Ab-initio Molecular Dynamics （VASP）

> [Ab-initio Molecular Dynamics (AIMD)](https://en.wikipedia.org/wiki/Molecular_dynamics#Ab-initio_molecular_dynamics) is a computational method that uses first principles to simulate the motion of atoms in a system. It is a type of molecular dynamics (MD) simulation that **does not rely on empirical potentials or force fields to describe the interactions between atoms, but rather calculates these interactions directly from the electronic structure of the system using quantum mechanics**.
> 
> In an ab initio MD simulation, the total energy of the system is calculated at each time step using density functional theory (DFT) or another method of quantum chemistry.  The forces acting on each atom are then determined from the gradient of the energy with respect to the atomic coordinates, and the equations of motion are solved to predict the trajectory of the atoms.



This document aims to record and summarize the ***computational details of AIMD simulation using VASP***.



## 1 Input

### 1.1 POSCAR

Preparing a  sufficiently large **supercell** (2$\times$2$\times$2, for example).

You can obtain it by following the description in: [Preparing a Super Cell - VASP Wiki](https://www.vasp.at/wiki/index.php/Preparing_a_Super_Cell), or using [VASPKIT](https://vaspkit.com) via the following command line:

```bash
echo -e "401\n1\n2 2 2\n"|vaspkit
mv SC222.vasp POSCAR
```

<div style="color:black; background-color:#7EC0EE; border: 1px solid #FFE0C3; border-radius: 10px; margin-bottom:0rem">
    <p style="margin:1rem; padding-left: 1rem; line-height: 2.5;">
        <b>💭Why do we use a supercell to perform MD simulations?</b><br/>

    <b>The size of the supercell imposes a limit on the maximum wavelength of lattice vibrations.</b> The supercell used in an MD simulation should be large enough to account for all vibration modes with significant contribution to the specific quantity of interest to be computed in MD. This can be estimated, e.g., from an appropriate phonon calculation, or from a series of MD simulations with different supercell sizes.<br/>

    Furthermore, in calculations considering for instance an adsorbate-substrate problem, or simulations of gases and liquids, the size of the unit cell should be large enough to remove unphysical interactions between atoms and their periodic images. Note that, the same holds also for relaxations of such systems.<br/>

    In summary, for your MD simulation, you should choose a supercell large enough to ensure an <a href="https://en.wikipedia.org/wiki/Ergodicity">ergodic simulation</a> and capture all long-wavelength vibrations of your system.<br/>
    </p>
</div>

> ***Note***:
> 
> If a continuation run is performed copy [CONTCAR](https://www.vasp.at/wiki/index.php/CONTCAR "CONTCAR") to [POSCAR](https://www.vasp.at/wiki/index.php/POSCAR "POSCAR") or possibly deliver initial velocities in the POSCAR file. They are written after the Wycoff positions in an own paragraph. If no initial velocities are provided random velocities are assumed at the beginning of the calculation. This is fully ok but the user should be aware that due to the initial random velocities the trajectories obtained from different calculations are difficult to compare.



### 1.2 INCAR

```
SYSTEM = Ab-initio MD

# ab initio
ISMEAR = 0        # Gaussian smearing
SIGMA  = 0.05     # smearing in eV

LREAL  = Auto     # projection operators in real space

ALGO   = Fast     # RMM-DIIS for electronic relaxation
PREC   = Normal   # precision
ISYM   = 0        # no symmetry imposed

# MD
IBRION = 0        # enable a molecular dynamics calculation
NSW    = 5000     # number of steps
POTIM  = 2.0      # MD time step in fs

MDALGO = 2        # Nosé-Hoover thermostat
ISIF   = 2
SMASS  = 1.0      # Nosé mass

TEBEG  = 1000     # temperature at beginning
TEEND  = 1000     # temperature at end

NCORE  = 2        # parallelization
```

[Thermostats](https://www.vasp.at/wiki/index.php/Category:Thermostats) are used in molecular-dynamics calculations within the [NVT ensemble](https://www.vasp.at/wiki/index.php/NVT_ensemble "NVT ensemble") and [NpT ensemble](https://www.vasp.at/wiki/index.php/NpT_ensemble "NpT ensemble") in order to apply a certain temperature to the ionic degrees of freedom.

<table cellpadding="5" cellspacing="0" border="1">
<tbody><tr>
<td></td>
<td colspan="4" style="text-align: center;"><a class="mw-selflink selflink">Thermostat</a>
</td></tr>
<tr>
<td style="text-align: center;"><a href="/wiki/index.php/Category:Ensembles" title="Category:Ensembles">Ensemble</a></td>
<td style="text-align: center;"><a href="/wiki/index.php/Andersen_thermostat" title="Andersen thermostat">Andersen</a></td>
<td style="text-align: center;"><a href="/wiki/index.php/Nose-Hoover_thermostat" title="Nose-Hoover thermostat">Nose-Hoover</a></td>
<td style="text-align: center;"><a href="/wiki/index.php/Langevin_thermostat" title="Langevin thermostat">Langevin</a></td>
<td style="text-align: center;"><a href="/wiki/index.php/MDALGO#MDALGO=13:_Multiple_Anderson_thermostats" title="MDALGO">Multiple Andersen</a>
</td></tr>
<tr>
<td style="text-align: center;"><a href="/wiki/index.php/NVE_ensemble" title="NVE ensemble">Microcanonical (NVE)</a></td>
<td colspan="4" style="text-align: center;"><a href="/wiki/index.php/MDALGO" title="MDALGO">MDALGO</a>=1,  <a href="/wiki/index.php/ANDERSEN_PROB" title="ANDERSEN PROB">ANDERSEN_PROB</a>=0.0
</td></tr>
<tr>
<td rowspan="2" style="text-align: center;"><a href="/wiki/index.php/NVT_ensemble" title="NVT ensemble">Canonical (NVT)</a></td>
<td style="text-align: center;"><a href="/wiki/index.php/MDALGO" title="MDALGO">MDALGO</a>=1</td>
<td style="text-align: center;"><a href="/wiki/index.php/MDALGO" title="MDALGO">MDALGO</a>=2</td>
<td style="text-align: center;"><a href="/wiki/index.php/MDALGO" title="MDALGO">MDALGO</a>=3</td>
<td style="text-align: center;"><a href="/wiki/index.php/MDALGO" title="MDALGO">MDALGO</a>=13
</td></tr>
<tr>
<td style="text-align: center;"><a href="/wiki/index.php/ISIF" title="ISIF">ISIF</a>=2</td>
<td style="text-align: center;"><a href="/wiki/index.php/ISIF" title="ISIF">ISIF</a>=2</td>
<td style="text-align: center;"><a href="/wiki/index.php/ISIF" title="ISIF">ISIF</a>=2</td>
<td style="text-align: center;"><a href="/wiki/index.php/ISIF" title="ISIF">ISIF</a>=2
</td></tr>
<tr>
<td rowspan="2" style="text-align: center;"><a href="/wiki/index.php/NpT_ensemble" title="NpT ensemble">Isobaric-isothermal (NpT)</a></td>
<td rowspan="2" style="text-align: center;">not available</td>
<td rowspan="2" style="text-align: center;">not available</td>
<td style="text-align: center;"><a href="/wiki/index.php/MDALGO" title="MDALGO">MDALGO</a>=3</td>
<td rowspan="2" style="text-align: center;">not available
</td></tr>
<tr>
<td style="text-align: center;"><a href="/wiki/index.php/ISIF" title="ISIF">ISIF</a>=3
</td></tr>
<tr>
<td style="text-align: center;"><a href="/wiki/index.php/NpH_ensemble" title="NpH ensemble">Isoenthalpic-isobaric (NpH)</a></td>
<td colspan="4" style="text-align: center;"><a href="/wiki/index.php/MDALGO" title="MDALGO">MDALGO</a>=3, <a href="/wiki/index.php/ISIF" title="ISIF">ISIF</a>=3, <a href="/wiki/index.php/LANGEVIN_GAMMA" title="LANGEVIN GAMMA">LANGEVIN_GAMMA</a>=<a href="/wiki/index.php/LANGEVIN_GAMMA_L" title="LANGEVIN GAMMA L">LANGEVIN_GAMMA_L</a>=0.0
</td></tr></tbody></table>

### 1.3 KPOINTS

```
K-Points
 0
Gamma
 1  1  1
 0  0  0
```

Since a sufficiently large super cell is used in MD, it is ok in this case to use only a single k-point in the calculations. Hence it is also possible to use the $\Gamma$-point only version which is significantly faster than the standard version.

### 1.4 POTCAR

*Pseudopotentials of the system*



## 2 Output

### 2.1 Pair correlation function

The pair correlation function is written out to the **[PCDAT](https://www.vasp.at/wiki/index.php/PCDAT "PCDAT")** file. The abscissa of that file is within mesh points of a selected grid and need to be converted to $\AA$. This is done by invoking the following short `awk` script on the command line:

```bash
awk <PCDAT >PCDAT.10ps ' NR==8 {pcskal=$1} NR==9 {pcfein=$1} NR>=13 {line=line+1; print (line-0.5)*pcfein/pcskal,$1} '
```

Also, you can do this via VASPKIT:

```bash
echo -e "725\n"|vaspkit
```

> Solids usually show sharp peaks in the pair correlation function since the ions only vibrate around fixed positions in the crystal lattice. In the liquid or amorphous state the distances are much more diffuse and one usually would expect no far order (but both can have near order).

### 2.2 Energy conservation

We can output the total energy for each molecular dynamics step by invoking the command:

```bash
grep "free  energy" OUTCAR|awk ' {print $5}' >> energy.dat
```







***Reference***:

- [Molecular dynamics calculations](https://www.vasp.at/wiki/index.php/Molecular_dynamics_calculations)

- [Molecular dynamics](https://www.vasp.at/tutorials/latest/md/)

- [Molecular dynamics - Tutorial](https://www.vasp.at/wiki/index.php/Molecular_dynamics_-_Tutorial)
