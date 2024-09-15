---
title: 'Linear Response (VASP)'
date: 2024-9-15
tags:
  - VASP
  - Linear Response
---

# Linear response

除了基态属性，VASP还可以计算系统如何对外部微扰作出反应。例如：

- 外部电场
- 原子位移
- 均匀应变

如果扰动很小，可以将其限制在一阶，则响应处于线性状态，可以通过**线性响应理论**进行计算。

线性响应的一个中心量是[**介电函数**](https://www.vasp.at/wiki/index.php/Category:Dielectric_properties#Dielectric_function)，它可以将外部电场与内部电位移联系起来。对原子位移的响应包括声子和电子—声子相互作用。利用均匀应变的响应，结合对外部电场的响应可以计算弹性张量和压电张量。

## 1 Polarization

周期性系统中的极化可以使用[极化的 berry 相公式](https://www.vasp.at/wiki/index.php/Berry_phases_and_finite_electric_fields#Modern_Theory_of_Polarization)（通常称为**现代极化理论**）计算。

> 严格而言，极化以及有限电场的应用都是基态属性，但是，因为它们可以用来计算静态介电张量、Born有效电荷和压电张量这些响应属性，因此在此处指出这种方法。

该方法通过设置 [`LCALCPOL`](https://www.vasp.at/wiki/index.php/LCALCPOL) 对系统施加一个[有限电场](https://www.vasp.at/wiki/index.php/Berry_phases_and_finite_electric_fields#Self-consistent_response_to_finite_electric_fields)计算电子极化。

## 2 [Static linear response](https://www.vasp.at/wiki/index.php/Static_linear_response:_theory)

微扰恒定的响应称为静态响应。不同的静态响应可以理解为**总能量相对于不同外部扰动的导数**。

> 对总能量$E$进行泰勒展开，得到
> $$
> \begin{aligned}
> E(u,\mathcal{E},\eta)=&E_{0}+ \\
> &\frac{\partial E}{\partial u_{m}}u_{m}+\frac{\partial E}{\partial\mathcal{E}_{\alpha}}\mathcal{E}_{\alpha}+\frac{\partial E}{\partial\eta_{j}}\eta_{j}+ \\
> &\frac{1}{2}\frac{\partial^{2}E}{\partial u_{m}\partial u_{n}}u_{m}u_{n}+\frac{1}{2}\frac{\partial^{2}E}{\partial\mathcal{E}_{\alpha}\partial\mathcal{E}_{\beta}}\mathcal{E}_{\alpha}\mathcal{E}_{\beta}+\frac{1}{2}\frac{\partial^{2}E}{\partial\eta_{j}\partial\eta_{k}}\eta_{j}\eta_{k}+ \\
> &\frac{\partial^{2}E}{\partial u_{m}\partial\mathcal{E}_{\alpha}}u_{m}\mathcal{E}_{\alpha}+\frac{\partial^{2}E}{\partial u_{m}\partial\eta_{j}}u_{m}\eta_{j}+\frac{\partial^{2}E}{\partial\mathcal{E}_{\alpha}\partial\eta_{j}}\mathcal{E}_{\alpha}\eta_{j}+\mathrm{terms~of~higher~order}
> \end{aligned}
> $$
>
> - 能量相对于**电场**的导数是**极化**：
>   $$
>   P_{\alpha}=-\frac{\partial E}{\partial\mathcal{E}_{\alpha}}\quad\mathrm{polarization}
>   $$
>
> - 能量相对于**原子位移**的导数是**力**：
>   $$
>   F_{m}=-\Omega_{0}\frac{\partial E}{\partial u_{m}}\quad\mathrm{forces}
>   $$
>
> - 能量相对于**晶格矢量**的导数是**应力张量**：
> $$
>   \sigma_{j}=\frac{\partial E}{\partial\eta_{j}}\quad\mathrm{stresses}
> $$
>   
> - 
>
> 由此引出***clamped-ion***或***frozen-ion***的定义：
> $$
> \begin{gathered}
> \bar{\chi}_{\alpha\beta}=-\frac{\partial^{2}E}{\partial\mathcal{E}_{\alpha}\partial\mathcal{E}_{\beta}}|_{u,\eta}\quad\mathrm{dielectric~susceptibility} \\
> \bar{C}_{jk}=\frac{\partial^{2}E}{\partial\eta_{j}\partial\eta_{k}}|_{u,\mathcal{E}}\quad\mathrm{elastic~tensor} \\
> \Phi_{mn}=\Omega_{0}\frac{\partial^{2}E}{\partial u_{m}\partial u_{n}}|_{\mathcal{E},\eta}\quad\mathrm{force-constants} \\
> \bar{e}_{\alpha k}=\frac{\partial^{2}E}{\partial\mathcal{E}_{\alpha}\partial\eta_{k}}|_{u}\quad\mathrm{piezoelectric~tensor} \\
> Z_{m\alpha}^{*}=-\Omega_{0}\frac{\partial^{2}E}{\partial u_{m}\partial\mathcal{E}_{\alpha}}|_{\eta}\quad\mathrm{Born~effective~charges} \\
> \Xi_{mj}=-\Omega_{0}\frac{\partial^{2}E}{\partial u_{m}\partial\eta_{j}}|_{\mathcal{E}}\quad\mathrm{force~response~internal~strain~tensor} 
> \end{gathered}
> $$
> 为了与实验结果进行比较，静态响应属性应考虑**离子弛豫**。通过以上泰勒展开并关注能量最小时的离子位置，得到
> $$
> \tilde{E}(\mathcal{E},\eta)=\min_{u}E(u,\mathcal{E},\eta)
> $$
> 物理***relaxed-ion***张量定义为：
> $$
> \begin{aligned}
> \chi_{\alpha\beta} &=\bar{\chi}_{\alpha\beta}+\Omega_{0}^{-1}Z_{m\alpha}^{*}(\Phi)_{mn}^{-1}Z_{n\beta}^{*} \quad\mathrm{dielectric~susceptibility} \\
> C_{jk} &=\bar{C}_{jk}+\Omega_{0}^{-1}\Xi_{mj}(\Phi)_{mn}^{-1}\Xi_{nk} \quad\mathrm{elastic~tensor} \\
> e_{\alpha j} &=\bar{e}_{\alpha j}+\Omega_{0}^{-1}Z_{m\alpha}^{*}(\Phi)_{mn}^{-1}\Xi_{nj} \quad\mathrm{piezoelectric~tensor} 
> \end{aligned}
> $$
> 每个方程右侧第二项称为**离子对介电极化率、弹性张量和压电张量的贡献**。
>
> - 离子对**介电张量**的贡献
>   $$
>   \epsilon_{ij}^{\mathrm{ion}}=\frac{4\pi}{\Omega}\sum_{kl}Z_{ik}^{*}\Phi_{kl}^{-1}Z_{lj}^{*}
>   $$
>
> - 离子对**弹性张量**的贡献
>   $$
>   C_{ik}^{\mathrm{ion}}=\sum_{kl}\Xi_{ij}\Phi_{jk}^{-1}\Xi_{kl}
>   $$
>
> - 离子对**压电张量**的贡献
>   $$
>   e_{ij}^{\mathrm{ion}}=\sum_{kl}Z_{ij}^{*}\Phi_{jk}^{-1}\Xi_{kl}
>   $$

### 2.1 Dielectric tensor

静态介电张量可以通过有限外电场极化的有限差分 [`LCALCEPS`](https://www.vasp.at/wiki/index.php/LCALCEPS) 或使用密度泛函微扰理论 [`LEPSILON`](https://www.vasp.at/wiki/index.php/LEPSILON) 计算。

> ⚠️***Note***: [`LEPSILON`](https://www.vasp.at/wiki/index.php/LEPSILON) 和 [`LCALCEPS`](https://www.vasp.at/wiki/index.php/LCALCEPS)  都能对介电张量产生相同的收敛结果，但是
>
> - **[`LEPSILON`](https://www.vasp.at/wiki/index.php/LEPSILON)只能用于局域或半局域交换关联泛函，并且适用于半导体和金属；**
> - **[`LCALCEPS`](https://www.vasp.at/wiki/index.php/LCALCEPS)可用于 [meta-GGA](https://www.vasp.at/wiki/index.php/Category:Meta-GGA) 或[杂化泛函](https://www.vasp.at/wiki/index.php/Category:Hybrid_functionals)，但[仅适用于具有带隙的体系](https://www.vasp.at/wiki/index.php/EFIELD_PEAD)。**

### 2.2 [Born effective charges - VASP Wiki](https://www.vasp.at/wiki/index.php/Born_effective_charges)

原子位移引起的极化变化在周期性系统中不是唯一定义的。周期性体系中的原子在不同的晶胞中重复，电荷可以广义化。Born有效电荷是定义这种**动态电荷**的一种方法。

Born有效电荷描述了体系的极化如何随离子的位置而变化，等价于一个电场诱导产生作用在离子上的力。一般而言，Born有效电荷是一个**张量**，也就是说，离子的极化方向和位移方向不一定平行。Born有效电荷有助于理解材料对外部刺激的响应，例如压电和铁电行为。

#### 2.2.1 Introduction

**动态电荷**定义为晶胞体积$\Omega_{0}$乘以宏观极化 ***P*** 在 *i* 方向上相对于原子$\kappa$亚晶格在 *j* 方向上的刚性位移的偏导数。

然而，极化在周期性系统中并不是唯一定义的，它依赖于由周期性边界条件固定的宏观电场$\mathcal{E}_{i}$。**Born有效电荷$Z^*$是宏观电场为零时极化对位置 *u* 的偏导数。**由于极化是总能量相对于宏观电场的一阶导数，所以$Z^{*}$可以用原子$\kappa$在 *j* 方向上的力 **F** 对$\mathcal{E}\_{i}$的偏导数来重新排列：
$$
Z_{\kappa,ij}^{*}=\frac{\Omega_{0}}{e}\frac{\partial\mathcal{P}_{i}}{\partial u_{\kappa,j}(q=0)}=\frac{1}{e}\frac{\partial F_{\kappa,j}}{\partial\mathcal{E}_{i}}\quad i,j=x,y,z
$$

> ***Mind***:
>
> - \*不表示复共轭，$Z^{*}$始终是实数；
> - $Z^*$在VASP中以 $\vert e\vert$ 为单位给出；
> - VASP输出$Z_{ij}^{*}$，其中 *i* 表示宏观电场，*j* 表示力的方向。在文献中，$Z_{ji}^{*}$很常见，即力方向 *j* 后跟电场方向 *i*。

#### 2.2.2 How to calculate

VASP实现了两种计算Born有效电荷的方法：

1. [`LCALCEPS = .TRUE.`](https://www.vasp.at/wiki/index.php/LCALCEPS) 

   通过沿三个笛卡尔方向施加[有限电场](https://www.vasp.at/wiki/index.php/Berry_phases_and_finite_electric_fields#Self-consistent_response_to_finite_electric_fields)并计算原子上的合力；

2. [`LEPSILON = .TRUE.`](https://www.vasp.at/wiki/index.php/LEPSILON)

   使用密度泛函微扰理论 （DFPT） 计算波函数相对于电场的导数。

以上方法可以与 [IBRION](https://www.vasp.at/wiki/index.php/IBRION) 结合使用，以获得额外的介电特性：

```
IBRION = 5 or 6 ! Calculated using finite differences.
IBRION = 7 or 8 ! Calculated using DFPT
```

Born 有效电荷将在 [OUTCAR](https://www.vasp.at/wiki/index.php/OUTCAR) 文件中给出：

```
BORN EFFECTIVE CHARGES (including local field effects) (in |e|, cummulative output)
BORN EFFECTIVE CHARGES (excluding local field effects) (in |e|, cummulative output)
```

> 局域场效应 - local field effects
>
> **局域场效应是指电场引起的轨道变化导致Hartree势和交换关联势的变化。**可以通过设置 [`LRPA`](https://www.vasp.at/wiki/index.php/LRPA) 以仅包含Hartree势的变化：
>
> ```
> LRPA = .TRUE.
> LCALCEPS = .TRUE. ! N.B. LEPSILON does not output the final Born effective charges.
> ```
>
> 这通常被称为**随机相位近似 (Random Phase Approximation, RPA) **中的响应，或“忽略局域场效应”。

### 2.3 Piezoelectric tensor

压电张量可以使用 [LCALCEPS](https://www.vasp.at/wiki/index.php/LCALCEPS) 通过[有限电场](https://www.vasp.at/wiki/index.php/Berry_phases_and_finite_electric_fields#Self-consistent_response_to_finite_electric_fields)的有限差分来计算，也可以使用 [LEPSILON](https://www.vasp.at/wiki/index.php/LEPSILON) 通过DFPT结合 `IBRION = 5,6 or 7,8` 计算。

### 2.4 Elastic tensor

通过 `IBRION = 5,6` 结合 `ISIF = 3` ，使用**应变有限差分**计算弹性张量。

> ⚠️***Note***: **由于在DFPT中没有实现应变扰动，因此不能使用 `IBRION = 7,8` 来计算弹性张量。**

### 2.5 Internal strain tensor

内部应变张量可以使用 `IBRION = 5,6` 通过有限差分计算，或者通过 `IBRION = 7,8` 使用DFPT计算。

## 3 Dynamic response

微扰随时间变化的响应称为动态响应。VASP中实现了[不同的方法](https://www.vasp.at/wiki/index.php/Category:Dielectric_properties#Dynamical_response:_Green-Kubo_and_many-body_perturbation_theory)和理论层次来计算***频率依赖的介电张量***。其中最简单的方法是使用***Green-kubo*公式**，通过 `LOPTICS` 设置。但是*Green-kubo*公式忽略了局域场效应，这意味着当频率为零（静态极限）时，如果`LRPA = .TRUE.`，只能通过DFPT或有限电场的有限差分重现计算。若要包含局域场效应，应使用 `ALGO = CHI`。此外，设置 `ALGO = BSE` 使用多体微扰理论框架中的Bethe-Salpeter方程计算介电函数，可以考虑电子-空穴相互作用。

# Calculation Details

提示：以下计算开始前均需要进行**结构优化**与**静态自洽**，并将静态自洽计算得到的[`WAVECAR`](https://www.vasp.at/wiki/index.php/WAVECAR)文件复制到相应的计算目录。

## 1 Static response with finite differences

`INCAR` 示例：

```
ISMEAR = 0
SIGMA  = 0.01
EDIFF  = 1E-6
ENCUT  = 400

! finite differences for polarization
LCALCPOL    = T 
LCALCEPS    = T   
LPEAD       = T     
EFIELD_PEAD = 0.01 0.01 0.01
IPEAD       = 4

! finite differences for ionic displacements
IBRION   = 6       
POTIM    = 0.015
NFREE    = 2

ISIF     = 3
```

> ⚠️**Notes**:
>
> - 电场太强可能会发生[齐纳隧穿 (Zener tunneling)](https://en.wikipedia.org/wiki/Zener_effect)，可以通过减少 $K$ 点数量、通过 [`EFIELD_PEAD`](https://www.vasp.at/wiki/index.php/EFIELD_PEAD) 设置减小电场；
> - [`LPEAD = .TRUE.`](https://www.vasp.at/wiki/index.php/LPEAD) 不支持金属体系。

## 2 Static dielectric response within DFPT

`INCAR` 示例：

```
SYSTEM = AlP

ISMEAR = 0
SIGMA  = 0.01
EDIFF  = 1E-6
ENCUT  = 400

LASPH  = T 

! finite differences for electric field
LPEAD       = T   

! dfpt for polarization
LEPSILON = T     

! dfpt for ionic displacements
IBRION   = 8
```

> ⚠️**Notes**:
>
> - 与 [LCALCEPS](https://www.vasp.at/wiki/index.php/LCALCEPS) 相比，[LEPSILON](https://www.vasp.at/wiki/index.php/LEPSILON) 适用于金属体系。
> - 有限差分和DFPT两种方法都没有对未占据状态求和。
> - 目前 [LEPSILON](https://www.vasp.at/wiki/index.php/LEPSILON) 不能用于明确依赖于轨道的交换关联泛函，例如类 HF 和杂化泛函。

## 3 Frequency-dependent dielectric response

完整的[频率相关介电函数](https://en.wikipedia.org/wiki/Dielectric_spectroscopy)包括**电子**和**离子**贡献。这些必须在 DFT 运行的基础上单独计算，然后相加：
$$
\varepsilon(\omega)=\varepsilon_{\mathrm{elec}}(\omega)+\varepsilon_{\mathrm{ion}}(\omega)
$$
**介电函数的虚部和实部通过[Kramers-Kronig 关系](https://en.wikipedia.org/wiki/Kramers–Kronig_relations)联系。**

从频率相关的介电函数出发，可以计算其他光学特性，例如反射率、吸收率和[光电导率](https://en.wikipedia.org/wiki/Optical_conductivity)。

Electron - `INCAR`:

```
ALGO   = Exact
ISMEAR = 0
SIGMA  = 0.01
EDIFF  = 1.E-6

LASPH  = T

! electronic dielectric function with Kramers-Kronig
LOPTICS  = T    ! frequency-dependent dielectric matrix
CSHIFT   = 0.1  ! complex shift (default)
NBANDS   = 32   ! sum over unoccupied bands
NEDOS    = 2000 ! frequency grid
```

Ion:dfpt - `INCAR`:

```
ALGO   = Normal
PREC   = High
ISMEAR = 0
SIGMA  = 0.01
EDIFF  = 1.E-8

! ionic dielectric function with Kramers-Kronig
! option 1: using purely DFPT
IBRION   = 8
LEPSILON = T
```

Ion:finite-differences - `INCAR`:

```
ALGO   = Normal
ISMEAR = 0
SIGMA  = 0.01
EDIFF  = 1.E-6

! ionic dielectric function with Kramers-Kronig
! option 2: using purely a finite differences approach
IBRION   = 6
NFREE    = 2
POTIM    = 0.015
TIME     = 0.1

LPEAD    = T
LCALCEPS = T
```

> ⚠️**Notes**:
>
> - 周期性体系的电子极化是$-fe/(2\pi)^3$乘以所有占据态的Berry相$\int_{BZ}\mathrm{d}\mathbf{k}\left\langle u_{n\mathbf{k}}^{(\mathcal{E}_{i})}\right|\mathrm{i}\nabla_{\mathbf{k}}\left|u_{n\mathbf{k}}^{(\mathcal{E}_{i})}\right\rangle $之和。
>
> - `LOPTICS = .TRUE.` 引入了未占据状态之和！因此必须设置 [`NBANDS`](https://www.vasp.at/wiki/index.php/NBANDS)以包含足够的未占据带。换言之，极化应该随着 `NBANDS` 的增加而收敛。
>
> - `POTCAR` 推荐采用[GW赝势](https://www.vasp.at/wiki/index.php/Available_pseudopotentials)。因为**在标准DFT计算中，未占据带对基态性质如总能量等没有贡献。为GW计算构建的赝势可以用于描述费米能级以上的性质。**因此，当增加 `NBANDS` 时，通常应该选择用于GW计算的赝势。
>
> - 声子频率与离子对介电函数的贡献有何关系？
>
>   如果一个声子模式是**dipole-active**的，它将出现在频率依赖的介电函数中。 如果声子模式中所有离子的电偶极矩在一个周期内大小发生变化，则这个声子模式是**dipole-active**的。这是动力学矩阵中特征向量的情况。
>
> - 单位转换：
>
>   频率的单位是
>   $$
>   [\omega]=\sqrt{\frac1{[m_{\mathrm{ion}}]}\left[\frac{\partial^2E}{\partial u_i\partial u_j}\right]}
>   $$
>   [`POMASS`](https://www.vasp.at/wiki/index.php/POMASS)以原子质量单位 (a.m.u.) 给出：
>   $$
>   [m_{\mathrm{ion}}]=1.660599\cdot10^{-27}\text{kg}
>   $$
>   将能量用 $eV$ 表示，离子位移用 $\mathring{A}$ 表示，则频率单位可以用 SI 表示为
>   $$
>   \begin{aligned}
>   1[\omega]&=\sqrt{\frac{\mathrm{eV/\mathring{A}}^2}{\mathrm{a.m.u.}}}\\&=\sqrt{\frac{1.602176487\cdot10^{-19}\mathrm{J/(10^{-20}m^{2})}}{1.660599\cdot10^{-27}\mathrm{kg}}}\\&=\sqrt{\frac{16.02176487\mathrm{J/m^{2}}}{1.660599\cdot10^{-27}\mathrm{kg}}}\\&=9.822517\cdot10^{13}\mathrm{s}^{-1}
>   \end{aligned}
>   $$
>   **其他典型频率单位**：
>
>   1. $eV$
>
>   $$
>   \begin{aligned}
>   x[E_{\mathrm{phonon}}]&=1[\omega]\cdot\hbar\\&=9.822517\cdot10^{13}s^{-1}\times1.054572\cdot10^{-34}\mathrm{Js}\\&=1.035855\cdot10^{-20}\mathrm{J}\\&=0.064652976 \mathrm{eV}
>   \end{aligned}
>   $$
>
>   2. $THz$
>
>   $$
>   \begin{aligned}
>   x[\nu]&=1[\omega]/2\pi\\&=98.22517/(2\pi)\cdot10^{12}\text{Hz}\\&=15.6330214\text{ THz}
>   \end{aligned}
>   $$
>
>   3. $cm^{-1}$
>
>   $$
>   \begin{aligned}
>   x[\lambda^{-1}]&=1[\nu]/c\\&=\frac{15.6330214\cdot10^{12}\mathrm{s}^{-1}}{29979245800\mathrm{~cm/s}}\\&=521.46146389\mathrm{~cm}^{-1}
>   \end{aligned}
>   $$
>
>   综合上述推导，得到
>
> $$
>   1\mathrm{~THz}=4.1357\mathrm{~meV}=33.356\mathrm{~cm}^{-1}\\1\mathrm{~meV}=0.242\mathrm{~THz}=8.066\mathrm{~cm}^{-1}\\1\mathrm{~cm}^{-1}=0.030\mathrm{~THz}=0.124\mathrm{~meV}.
> $$
>







***Ref***:

- [Linear response - VASP Wiki](https://www.vasp.at/wiki/index.php/Category:Linear_response)
- [Static linear response: theory - VASP Wiki](https://www.vasp.at/wiki/index.php/Static_linear_response:_theory)
- [*Phys. Rev. B* **58**, 6224 (1998)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.58.6224)
- [*Phys. Rev. B* **55**, 10355 (1997)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.55.10355)
- [*Phys. Rev. B* **72**, 035105 (2005)](https://doi.org/10.1103/PhysRevB.72.035105)
- [Born effective charges - VASP Wiki](https://www.vasp.at/wiki/index.php/Born_effective_charges)
- [*Phys. Rev. B* **73**, 045112 (2006)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.73.045112)
- [Linear response - Part 1](https://www.vasp.at/tutorials/latest/response/part1/)
