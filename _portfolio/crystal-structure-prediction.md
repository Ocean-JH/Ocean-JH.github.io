---
title: "Crystal Structure Prediction"
excerpt: "Explore the configuration space of Ge-Sb-Te system using genetic algorithm , aiming to find new (meta)stable materials.<br/><img src='/images/potential-energy-surface.png'>"
collection: portfolio
---

# 1. CSP with Genetic Algorithms

Utilizing a combination of genetic algorithms and first-principles methods, I conducted crystal structure prediction studies on the Ge-Sb-Te system, a pivotal phase change material system, to uncover its inherent complexities and expand our knowledge of its configuration space.

In this work, ***crossover*** and a variety of ***mutation*** operators are used to improve the sampling efficiency of the algorithm, including *mutation*, *permutation*, *transmutation*, etc. The energy above hull($E_{hull}$) is used as **Fitness function** for population iteration. After 50 generations of iteration, the resulting candidate structures are sorted by fitness, and then further validated according to the following carefully considered screening criteria:

1. Structures with formation energy($E_{form}$) < 0 and energy above hull($E_{hull}$) < 0.1 eV/atom is selected to calculate the *elastic constant*, and the stability of the structure is further determined by [***Born elastic criterion***](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.224104) (*Phys. Rev. B* **90**, 224104 (2014)) and frequency;
2. Structures is preliminary deduplicated by the number of atoms of each element and space group;
3. Calculate the *phonon spectrum /AIMD* of the structures to determine the dynamic stability of the material;
4. After the above screening steps, *Born effective charge*, *density of state*, *band structure*, *dielectric function* and *transport properties* of the remaining structure are further analyzed to identify potential application scenarios.

The objective of this work is to enhance comprehension of this intriguing system and potentially identify novel structures or compositions conducive to phase-change storage applications.

# 2. CSP with structural prototype

Given the extensive *crystal material database*, we can search for new materials through ***element substitution*** based on ***structural prototype***. I adopted this strategy for the impressive ***GeTe-Sb2Te3 pseudo-binary line*** of the Ge-Sb-Te system, due to the large number of phase-change storage materials found along it.
The overall process is similar to the above work, with a more stringent energy screening criterion ($E_{form} < E_{cubic}$) due to the prevalence of ***metastable cubic*** phases in structures along pseudo-binary line. In addition, since there are a significant number of structures in the database that share the same prototype, structure de-duplication is a priority.