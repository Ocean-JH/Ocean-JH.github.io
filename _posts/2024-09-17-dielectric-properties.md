---
title: 'Dielectric Properties (VASP)'
date: 2024-9-17
tags:
  - VASP
  - Dielectric function
---

# Dielectric properties

## 1 Dielectric function

当外部电场 ***E*** 作用于介质时，电子和离子电荷都会对扰动场产生反应。对于介电材料，一种简单的方法是认为束缚电荷会在介质内部产生偶极子，从而产生感应极化 ***P***。两个场的综合作用表示为电位移场 ***D***，由下式给出：

$$
D=E+4\pi P
$$

**如果外场强度不足以显著改变介电介质的性质，则可以在所谓的线性响应范围内处理感应极化。**此处介电介质对外部场的反应信息由**介电函数**给出：

$$
\epsilon_{\alpha\beta}=\delta_{\alpha\beta}+4\pi\frac{\partial P_{i}}{\partial E_{j}}
$$

假设体系具有时间反演对称性，这将导致

$$
D_\alpha(\omega)=\epsilon_{\alpha\beta}(\omega)E_\beta(\omega)
$$

根据外场的性质，有不同的方法计算 $\epsilon$：

- 如果 ***E*** 是静态的，可以采用基于有限差分的微扰方法或者密度泛函微扰理论（DFPT）。
- 如果扰动是一个时间依赖的 ***E***，例如在测量光的吸收率、反射率、[磁光克尔效应（Magneto-optical Kerr effect, MOKE）](https://en.wikipedia.org/wiki/Magneto-optic_Kerr_effect)等时，响应将取决于外场的频率。对于这些情况，必须采用基于时变线性响应（例如，***Green-Kubo***）或多体微扰理论的方法。

### 1.1 Static response

#### [LEPSILON](https://www.vasp.at/wiki/index.php/LEPSILON): density-functional-perturbation theory (DFPT)

设置 `LEPSILON=.True.`，VASP使用DFPT计算包含或不包含局域场效应（[`LRPA`](https://www.vasp.at/wiki/index.php/LRPA)）的静态 *ion-clamped* 介电张量。导数的计算使用Sternheimer方程，避免波函数周期性部分导数的显式计算。此方法不需要通过 [`NBANDS`](https://www.vasp.at/wiki/index.php/NBANDS) 包含空状态。

计算结束后，包含（`LRPA = .True.`）或不包含（`LRPA = .False.`）局域场效应的 $\epsilon$ 值输出在 `OUTCAR` 文件中。通过比较不包含局域场效应的值和 $\epsilon$ 的静态极限$\lim_{\omega\to0}\epsilon(\omega)$（通过 `LOPTICS = .True.` 获得）执行一致性检查。

#### [LCALCEPS](https://www.vasp.at/wiki/index.php/LCALCEPS): finite differences approach

设置 `LCALCEPS = .True.`，介电张量根据极化的导数计算得到：

$$
\epsilon_{ij}^{\infty}=\delta_{ij}+\frac{4\pi}{\epsilon_{0}}\frac{\partial P_{i}}{\partial\mathcal{E}_{j}}\quad i,j=x,y,z.
$$

此处的导数通过有限差分显式计算得到。扰动电场的方向和强度必须在 `INCAR` 文件中使用 [`EFIELD_PEAD`](https://www.vasp.at/wiki/index.php/EFIELD_PEAD) 标签指定。与DFPT一样，计算结束后VASP将在 `OUTCAR` 文件中写入介电张量。考虑局域场效应的控制通过变量 [`LRPA`](https://www.vasp.at/wiki/index.php/LRPA) 完成。

### 1.2 Dynamical response: Green-Kubo and many-body perturbation theory

#### [LOPTICS](https://www.vasp.at/wiki/index.php/LOPTICS): Green-Kubo formula

计算得到电子基态后，可以通过 `LOPTICS` 评估频率依赖的介电函数。$\epsilon$ 的**虚部采用显式方程通过对空状态求和确定**：

$$
\epsilon_{\alpha\beta}^{(2)}\left(\omega\right)=\frac{4\pi^{2}e^{2}}{\Omega}\mathrm{lim}_{q\to0}\frac{1}{q^{2}}\sum_{c,v,\mathbf{k}}2w_{\mathbf{k}}\delta(\epsilon_{c\mathbf{k}}-\epsilon_{v\mathbf{k}}-\omega)\times\langle u_{c\mathbf{k}+\mathbf{e}_{\alpha}q}|u_{v\mathbf{k}}\rangle\langle u_{v\mathbf{k}}|u_{c\mathbf{k}+\mathbf{e}_{\beta}q}\rangle
$$

此处的索引 *c* 和 *v* 分别指导带和价带状态，$u_{c\mathbf k}$ 是轨道在 *k-point* ***k***处的晶胞周期性部分。

**实部通过Kramers-Kronig变换求解**：

$$
\epsilon_{\alpha\beta}^{(1)}(\omega)=1+\frac{2}{\pi}P\int_{0}^{\infty}\frac{\epsilon_{\alpha\beta}^{(2)}(\omega^{\prime})\omega^{\prime}}{\omega^{\prime2}-\omega^{2}+i\eta}d\omega^{\prime}
$$

复位移 $\eta$ 由参数 [`CSHIFT`](https://www.vasp.at/wiki/index.php/CSHIFT) 确定。

在这一近似中忽略了局域场效应，即势的晶胞周期性部分的变化。该效应可以通过密度泛函微扰理论(`LEPSILON = .True.`)或GW计算进行评估。

该方法需要两步：

1. 首先，获得电子基态；
2. 然后，增加 `INCAR` 文件中 `NBANDS` 的值，以包含未占据状态。

此外，`INCAR` 还应包括 [`CSHIFT`](https://www.vasp.at/wiki/index.php/CSHIFT) 值(应用于洛伦兹函数代替$\delta$-函数的展宽)和 [`NEDOS`](https://www.vasp.at/wiki/index.php/NEDOS)（$\omega$的频率网格点数量）

> **Notes**:
>
> - 时刻检查涉及未占据状态数的收敛性。
> - 使用 `LOPTICS = .TRUE.` 的方法需要相当数量的空导带状态。**通常只有当 [`INCAR`](https://www.vasp.at/wiki/index.php/INCAR) 文件中的参数 [`NBANDS`](https://www.vasp.at/wiki/index.php/NBANDS) 相对于VASP默认值大约两倍或三倍时，才能获得合理的结果。 **
> - 在许多情况下，最好显著增加 [`NEDOS`](https://www.vasp.at/wiki/index.php/NEDOS) 值。**强烈建议使用 `NEDOS = 2000` 左右的值。**



🔖[LEPSILON](https://www.vasp.at/wiki/index.php/LEPSILON) *Vs* [LOPTICS](https://www.vasp.at/wiki/index.php/LOPTICS) - Pros & Cons:

- Pros of [LEPSILON](https://www.vasp.at/wiki/index.php/LEPSILON):
  - 不需要导带
  - [在RPA和DFT层面包含局域场效应](https://www.vasp.at/wiki/index.php/ACFDT/RPA_calculations)
- Cons of [LEPSILON](https://www.vasp.at/wiki/index.php/LEPSILON):
  - 目前只适用于静态属性
  - 需要相对耗时的迭代过程
  - 不支持HF或杂化泛函

不建议在单次计算中同时设置 `LOPTICS = .True.` 和 `LEPSILON = .True.`。密度泛函微扰理论（`LEPSILON = .True.`）不需要增加 [NBANDS](https://www.vasp.at/wiki/index.php/NBANDS)，否则速度会慢得多。

#### [ALGO](https://www.vasp.at/wiki/index.php/ALGO) = TDHF: Casida equation

该选项执行[时间依赖Hartree-Fock](https://www.vasp.at/wiki/index.php/Bethe-Salpeter-equations_calculations#Time-dependent_Hartree-Fock_calculation)或[时间依赖的密度泛函理论（TDDFT）](https://www.vasp.at/wiki/index.php/Bethe-Salpeter-equations_calculations#Time-dependent_DFT_calculation)计算。它遵循Casida方程并使用时间演化偶极子的傅里叶变换来计算 $\epsilon$。

[`NBANDS`](https://www.vasp.at/wiki/index.php/NBANDS) 的数量控制时间演化中出现带的数量。与 [`LOPTICS`](https://www.vasp.at/wiki/index.php/LOPTICS) 相比该方法可以减少空状态。

时间依赖的核心选择由 [`AEXX`](https://www.vasp.at/wiki/index.php/AEXX)，[`HFSCREEN`](https://www.vasp.at/wiki/index.php/HFSCREEN) 和 [`LFXC`](https://www.vasp.at/wiki/index.php/LFXC) 标签控制。

- 对于使用杂化泛函的计算，[AEXX](https://www.vasp.at/wiki/index.php/AEXX) 控制交换关联势中直接交换的比例，[`HFSCREEN`](https://www.vasp.at/wiki/index.php/HFSCREEN) 指定 [range-separated hybrid functionals](https://www.vasp.at/wiki/index.php/Hybrid_functionals:_formalism) 中的range-separation参数。

- 对于纯TDDFT计算，[`LFXC`](https://www.vasp.at/wiki/index.php/LFXC) 在时间演化方程中使用局域交换关联核：
  
  $$
  f_{\mathrm{xc}}^{\mathrm{loc}}\left(\mathbf{r},\mathbf{r}'\right)=\frac{\delta^{2}\left\{E_{\mathrm{c}}^{\mathrm{DFT}}+\left(1-c_{\mathrm{x}}\right)E_{\mathrm{x}}^{\mathrm{DFT}}\right\}}{\delta\rho(\mathbf{r})\delta\rho\left(\mathbf{r}'\right)}
  $$
  
  其中 $c_{\mathrm x}$ 是由 [`AEXX`](https://www.vasp.at/wiki/index.php/AEXX) 设置的交换相互作用比例。

#### [ALGO](https://www.vasp.at/wiki/index.php/ALGO) = TIMEEV: delta-pulse electric field

利用delta脉冲电场探测所有跃迁，并根据偶极动量的[时间演化](https://www.vasp.at/wiki/index.php/Time_Evolution)计算介电函数。通过 [`LHARTREE = .True.`](https://www.vasp.at/wiki/index.php/LHARTREE) 和 [`LFXC = .True.`](https://www.vasp.at/wiki/index.php/LFXC) 设置正确的时间依赖核，该算法能够完全再现标准[Bethe-Salpeter 计算](https://www.vasp.at/wiki/index.php/Bethe-Salpeter_equations)得到的吸收光谱。

时间步长由 [`CSHIFT`](https://www.vasp.at/wiki/index.php/CSHIFT) 和 [`PREC`](https://www.vasp.at/wiki/index.php/PREC) 自动控制。这意味着若 `CSHIFT` 值越小，精度选择越精确，则VASP执行的时间步数就越多，计算的成本也就越高。

参与时间传播的价带和导带数量分别由 [`NBANDSO`](https://www.vasp.at/wiki/index.php/NBANDSO) 和 [`NBANDSV`](https://www.vasp.at/wiki/index.php/NBANDSV) 标签设置。选择带隙附近较少数量的带来重现光学测量。

最后，在傅里叶变换和计算频率依赖介电函数时使用的最大能量由 [`OMEGAMAX`](https://www.vasp.at/wiki/index.php/OMEGAMAX) 设定，频率网格的采样由 [`NEDOS`](https://www.vasp.at/wiki/index.php/NEDOS) 控制。

#### [ALGO](https://www.vasp.at/wiki/index.php/ALGO) = CHI: polarizability within RPA approximation

通过随机相位近似（Random-Phase approximation, RPA）计算频率介电函数。VASP通过在 `INCAR` 文件设置 `ALGO = CHI` 计算极化率 $\chi$，然后使用下式计算介电函数：

$$
\epsilon_{\mathbf{GG}'}^{-1}(\mathbf{q},\omega)=\delta_{\mathbf{GG}'}+v(\mathbf{q}+\mathbf{G})\chi_{\mathbf{GG}'}(\mathbf{q},\omega)
$$

此处 *v* 是描述电子—电子相互作用的基本库伦势。这需要增加 `NBANDS` 以包含未占据态。

计算极化率有两种办法：

1. 设置 [LSPECTRAL = .True.](https://www.vasp.at/wiki/index.php/LSPECTRAL)，VASP将避免直接计算 $\chi$，而是使用一个快速的矩阵—向量乘积。然而对于某些 `CSHIFT` 值和 `NOMEGA` 值，这样做可能会在低频处引入伪峰。
2. 设置 [LSPECTRAL = .False.](https://www.vasp.at/wiki/index.php/LSPECTRAL)，直接计算 $\chi$，但比前一种方法要慢得多。

#### [ALGO](https://www.vasp.at/wiki/index.php/ALGO) = BSE: macroscopic dielectric function including excitons

通过求解[Bethe-Salpeter方程](https://www.vasp.at/wiki/index.php/Category:Bethe-Salpeter_equations)计算宏观介电函数 $\epsilon_M$。将电子-空穴对视为一种新的准粒子，称为***激子***，并利用特征向量 $X_\lambda^{cv\mathbf{k}}$ 和特征值 $\omega_\lambda$ 构建介电函数：

$$
\epsilon_M(\mathbf{q},\omega)=1+v(\mathbf{q})\sum_{\lambda\lambda^{\prime}}\sum_{c,v,\mathbf{k}}\sum_{c^{\prime},v^{\prime},\mathbf{k}^{\prime}}\langle c\mathbf{k}|e^{i\mathbf{q}\mathbf{r}}|v\mathbf{k}\rangle X_\lambda^{cv\mathbf{k}}\langle c^{\prime}\mathbf{k}^{\prime}|e^{-i\boldsymbol{q}\mathbf{r}}|v^{\prime}\mathbf{k}^{\prime}\rangle X_{\lambda^{\prime}}^{c^{\prime}v^{\prime}\mathbf{k}^{\prime},*}\times S_{\lambda,\lambda^{\prime}}^{-1}\left(\frac{1}{\omega_\lambda-\omega-i\delta}+\frac{1}{\omega_\lambda+\omega+i\delta}\right)
$$

包含在BSE哈密顿量中已占据和未占据状态的数量分别由 [NBANDSO](https://www.vasp.at/wiki/index.php/NBANDSO) 和 [NBANDSV](https://www.vasp.at/wiki/index.php/NBANDSV) 控制。通常情况下，只需要带隙上下的几个带就可以收敛光谱，且**随着带数量的增加，对内存的要求也会迅速增加**。

在与光学实验（如吸收、反射、MOKE）的比较方面，***q*** 是光子动量。通常考虑 $q\to0$ 的极限。此外，在所谓[Tamm-Dancoff近似](https://www.vasp.at/wiki/index.php/Bethe-Salpeter_equations#Theory#Tamm-Dancoff_approximation)中，谐振项和反谐振项之间的耦合可以被关闭。该近似可以通过设置变量 [`ANTIRES = 0`](https://www.vasp.at/wiki/index.php/ANTIRES) 激活，设置 [`ANTIRES = 1 or 2`](https://www.vasp.at/wiki/index.php/ANTIRES) 将考虑耦合，但会增加计算成本。

## 2 Level of approximation

### 2.1 Microscopic and macroscopic quantities

实验上测量一种性质时，实验数据呈现的是宏观量；而另一方面，计算上更容易获取微观量。**为了比较实验和计算结果，微观量（例如介电函数）必须在多次重复的单胞上取平均。**通过下式体现宏观介电函数 $\epsilon_M(\mathbf{q},\omega)$ 和微观介电函数的关系：

$$
\epsilon_M(\mathbf{q},\omega)=\frac{1}{\epsilon_{\mathbf{G}=0,\mathbf{G}^{\prime}=0}^{-1}(\mathbf{q},\omega)}
$$

其中 $\epsilon_{\mathbf{G}=0,\mathbf{G}^{\prime}=0}^{-1}(\mathbf{q},\omega)$ 是 $\mathbf{G}=\mathbf{G}^{\prime}=0$ 处的逆介电函数。

> **Note**：
>
> 这并不意味着 $\epsilon_{M}(\mathbf{q},\omega)=\epsilon_{\mathbf{G}=0,\mathbf{G}^{\prime}=0}(\mathbf{q},\omega)$。整个矩阵 $\epsilon_{M}(\mathbf{q},\omega)$ 必须求逆且它是在 $\mathbf{G}=\mathbf{G}^{\prime}=0$ 处用于计算 $\epsilon_M(\mathbf{q},\omega)$ 的分量。

### 2.2 Finite momentum dielectric function

在光学极限下，入射光子 ***q*** 的动量几乎为0，因为电场的波长比单胞的尺寸大好几倍。由于库伦势在非常小的动量下发散，介电函数的光学极限必须通过 $q\to0$ 的极限得到。例如，在完全BSE的独立粒子近似情况下，得到

$$
\lim_{\mathbf{q}\to0}\frac{\langle c\mathbf{k}+\mathbf{q}|e^{\mathrm{i}\mathbf{q}\cdot\mathbf{r}}|v\mathbf{k}\rangle}{q}\approx\lim_{\mathbf{q}\to0}\frac{\langle c\mathbf{k}+\mathbf{q}|1+\mathrm{i}\mathbf{q}\cdot\mathbf{r}|v\mathbf{k}\rangle}{q}=\hat{\mathbf{q}}\cdot\langle c\mathbf{k}+\mathbf{q}|\mathbf{r}|v\mathbf{k}\rangle
$$

VASP还可以分析有限动量激子的影响。为计算有限动量下的吸收光谱，将 [KPOINT_BSE](https://www.vasp.at/wiki/index.php/KPOINT_BSE) 设置为所需 ***q*** 点的索引。

### 2.3 Local fields in the Hamiltonian

局域场，即具有有限 ***G*** 的项，在评估极化率时可以在库仑势中打开或关闭。VASP将在 `OUTCAR` 文件中区分两种情况的结果：

```
MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects)
```

```
BORN EFFECTIVE CHARGES (including local field effects)
```

```
PIEZOELECTRIC TENSOR (including local field effects)
```

和

```
MACROSCOPIC STATIC DIELECTRIC TENSOR (excluding local field effects)
```

```
BORN EFFECTIVE CHARGES (excluding local field effects)
```

```
PIEZOELECTRIC TENSOR (excluding local field effects)
```

还可以采用另一种近似以在评估极化率时忽略交换关联核的贡献。这等价于所谓的随机相位近似（random-phase approximation, RPA），可以通过在 `INCAR` 中设置 [`LRPA = .True.`](https://www.vasp.at/wiki/index.php/LRPA) 激活。

### 2.4 Ion-clamped vs relaxed-ion/dressed dielectric function

在静态或动态响应状态下计算的介电函数都没有考虑入射电场引起离子位置变化而产生的影响。这可以通过计算 *relaxed-ion* （或 *dressed*）介电函数 $\bar{\epsilon}$ 进行修正：

$$
\bar{\epsilon}_{\alpha\beta}=\epsilon_{\alpha\beta}+\Omega_{0}^{-1}Z_{m\alpha}(\Phi^{-1})_{mn}Z_{n\beta},
$$

其中 $\Omega_{0}$ 是单胞的体积，$Z$ 是Born有效电荷，$\Phi_{mn}$ 是力常数矩阵。

### 2.5 Density-density versus current-current response functions

哈密顿量中是否包含电磁场取决于规范的选择。例如，经典电场既可以用一个标量势 $\phi$ 描述，也可以用纵向矢量势 ***A*** 在不完全Weyl规范（$\phi=0$）下描述。前者意味着扰动势对电子密度的耦合，而后者意味着矢量势对电流的耦合。这样的基本结果是我们可以定义两种不同的响应函数：

- 密度—密度响应函数 - $\chi_{\rho\rho}$
- 电流—电流响应函数 - $\chi_{jj}$

在比较周期性系统的实验和计算光学性质时，这种情况是一个常见的误差来源。

事实上，与纵向场相关的微扰由密度—密度极化率函数 $\chi_{\rho\rho}$ 描述，例如经典极限下的激光场；而横向场由电流—电流极化率函数 $\chi_{jj}$ 描述，它实际上是一个 $3\times3$ 张量，例如获取MOKE所需的张量。基本上，时间依赖的密度只通过连续性方程与电流的纵向部分相关联。可以通过下式连接这两种响应函数：

$$
q^2\chi_{jj}(\mathbf{q},\omega)=\omega^2\chi_{\rho\rho}(\mathbf{q},\omega)
$$

这保证了通过两种方法得到的介电函数在有限动量和频率下相匹配：$\epsilon[\chi_{\rho\rho}]=\epsilon[\chi_{jj}]$。

电流—电流介电函数在 $q\to0$ 和 $q=0$ 处都是准确的。特别是对于金属，它将重现 $\omega=0$ 处Drude tail的正确行为。然而，$\chi_{jj}$ 更易出现数值不稳定。

## 3 Other dielectric properties

频率依赖的线性光谱，例如折射率 $n(\omega)$，消光系数 $k(\omega)$，吸收系数 $\alpha(\omega)$，能量损失函数 $L(\omega)$ 和反射率 $R(\omega)$ 可以通过介电函数的实部 $\varepsilon_{1}(\omega)$ 和虚部 $\varepsilon_{2}(\omega)$ 计算：

$$
\begin{aligned}
n(\omega)&=\left(\frac{\sqrt{\varepsilon_{1}^{2}+\varepsilon_{2}^{2}}+\varepsilon_{1}}{2}\right)^{\frac{1}{2}}\\
k(\omega)&=\left(\frac{\sqrt{\varepsilon_1^2+\varepsilon_2^2}-\varepsilon_1}2\right)^{\frac12}\\
\alpha(\omega)&=\frac{\sqrt2\omega}{c}\left(\sqrt{\varepsilon_1^2+\varepsilon_2^2}-\varepsilon_1\right)^{\frac12}\\
L(\omega)&=\mathrm{Im}\left(\frac{-1}{\varepsilon(\omega)}\right)=\frac{\varepsilon_{2}}{\varepsilon_{1}^{2}+\varepsilon_{2}^{2}}\\
R(\omega)&=\frac{(n-1)^2+k^2}{(n+1)^2+k^2}
\end{aligned}
$$

### 3.1 Electron energy loss spectroscopy (EELS)

在EELS实验中，一束具有明确能量的窄电子束射向样品，然后这些电子通过激发等离激元、电子—空穴对或其他高阶准粒子而损失能量给样品。损失函数可以表示为：

$$
\mathrm{EELS}=-\mathrm{Im}\left[\epsilon^{-1}(\omega)\right]
$$

### 3.2 Optical conductivity

从麦克斯韦方程组和欧姆定律的微观形式可以得出张量介电函数和光电导率 $\sigma(\omega)$ 之间的关系式：

$$
\sigma_{\alpha\beta}(\omega)=\mathrm{i}\frac{\omega}{4\pi}\left[\delta_{\alpha\beta}-\epsilon_{\alpha\beta}(\omega)\right]
$$

### 3.3 Optical absorption

对于穿过介质的电磁波，可以将电场表示为 $\mathbf{E}(\mathbf{r},t)=\mathbf{E}_{0}e^{-\mathrm{i}(\omega t-\mathbf{q}\cdot r)}$，而介质对波传播的影响包含在色散关系 $\omega=\omega(\mathbf{q})$ 中。使用麦克斯韦方程组，得到

$$
q^2=\frac{\omega^2}{c^2}\epsilon(\omega)
$$

假设材料的磁导率等于真空磁导率，由上式可知，折射率可以表示为 $n=\sqrt{\epsilon(\omega)}=\tilde{n}+\mathrm{i}k$。由于 $n$ 是复数，$\mathbf{E}(\mathbf{r},t)$ 的指数因子会有一个阻尼系数 $e^{-\frac{\omega}{c}k\hat{q}\cdot\mathbf{r}}$，这个阻尼系数解释了介质对电磁能量的吸收。根据这个关系可以定义吸收系数：

$$
\alpha(\omega)=\frac{2\omega}{c}k(\omega)
$$

### 3.4 Reflectance

承上一小节，可以定义法向入射时的反射率系数：

$$
R=\frac{(1-\tilde{n})^2+k^2}{(1+\tilde{n})^2+k^2}
$$

该方程可以推广到任意入射角 $\theta$，从而得到Fresnel方程的一般形式。

### 3.5 Magneto-optical Kerr effect (MOKE)

磁光克尔效应是指入射电磁波与材料的有限磁矩相互作用。相互作用一般与介质的磁化有关，但也有可以观察到有限MOKE的反铁磁系统。通常，反射波相对于入射 ***E**-field* 会获得一个额外的复相位。对于表面或二维材料（如六方$BN$、$MoS_2$），该相位可以使用电流—电流介电张量的非对角分量来计算：

$$
\theta_\mathrm{K}(\omega)=-\mathrm{Re}\left[\frac{\epsilon_{xy}(\omega)}{(\epsilon_{xx}(\omega)-1)\sqrt{\epsilon_{xx}(\omega)}}\right]
$$

## 4 Electric response combined with perturbations of the ionic degrees of freedom

### 4.1 Low-frequency corrections from atomic displacements

离子运动到低频区域的修正可以通过下式添加到 $\epsilon_{\alpha\beta}^{\infty}$：

$$
\epsilon_{\alpha\beta}(\omega)=\epsilon_{\alpha\beta}^{\infty}+\frac{4\pi e^{2}}{\Omega_{0}}\sum_{\nu}\frac{S_{\alpha\beta,\nu}}{\omega_{\nu}^{2}-(\omega+\mathrm{i}\eta)^{2}}
$$

其中 $\omega_{\nu}$ 是模式 $\nu$ 的声子频率，$S_{\alpha\beta,\nu}$ 是由下式定义的mode-oscillator强度：

$$
S_{\alpha\beta,\nu}=\left(\sum_{I,\delta}Z_{I\alpha\delta}^*\varepsilon_{I\delta,\nu}^*(\mathbf{q}=\mathbf{0})\right)\left(\sum_{J,\delta^{\prime}}Z_{J\beta\delta^{\prime}}^*\varepsilon_{J\delta^{\prime},\nu}(\mathbf{q}=\mathbf{0})\right)
$$

此处 $Z_{J\beta\delta^{\prime}}^*$ 是Born有效电荷，而 $\varepsilon_{J\delta^{\prime},\nu}(\mathbf{q}=\mathbf{0})$ 是原子 *J* 沿方向 $\delta^{\prime}$ 与振动模式 $\nu$ 相关的本征位移。

> 关于**声子频率**和**特征位移**计算的理论和方法，详见 [Phonons](https://www.vasp.at/wiki/index.php/Phonons:_Theory) 页面。

在 `INCAR` 文件中，可以通过 `IBRION = 5, 6` （有限差分）或 `IBRION = 7, 8` （DFPT）并设置 [`LEPSILON = .True.`](https://www.vasp.at/wiki/index.php/LEPSILON) 或 [`LCALCEPS = .True.`](https://www.vasp.at/wiki/index.php/LCALCEPS) 以考虑低频校正。

> ***Polar materials***:
>
> 对于极性材料，在 $\Gamma$ 附近存在不连续性，即 $\omega_{\nu}^{2}(\mathbf{q}\to\mathbf{0})\neq\omega_{\nu}^{2}(\mathbf{q}=\mathbf{0})$。事实上，对于一个给定的酉方向矢量 ***q***，可以证明Lyddane-Sachs-Teller关系成立：
> 
> $$
> \prod_{\nu}\frac{\omega_{\nu}^{2}(\mathbf{q}\to\mathbf{0})-\omega^{2}}{\omega_{\nu}^{2}(\mathbf{q}=\mathbf{0})-\omega^{2}}=\frac{\sum_{\alpha\beta}q_{\alpha}\epsilon_{\alpha\beta}(\omega)q_{\beta}}{\sum_{\alpha\beta}q_{\alpha}\epsilon_{\alpha\beta}^{\infty}q_{\beta}}
> $$
> 
> 这意味着零动量下LO和TO模式之间的频率劈裂延续了对介电函数的评估。
>
> 为获得平滑的声子色散，并在评估介电函数的光学极限时正确考虑LO-TO劈裂，详见 [LO-TO splitting](https://www.vasp.at/wiki/index.php/Phonons:_Theory#Long-range_interatomic_force_constants_(LO-TO_splitting)) 页面。

### 4.2 Corrections from strain

介电张量可以包含在弹性张量 $C_{jk}$ 的评估中（详见[静态线性响应理论](https://www.vasp.at/wiki/index.php/Static_linear_response:_theory)）。虽然这个量通常在固定 ***E**-field* 下评估，但在绝缘材料层间放置薄膜的情况下，计算固定位移场 ***D*** 的弹性张量更加方便，因为边界条件将该矢量的分量固定在与表面垂直的方向上。假设 $C_{jk}^E$ 是固定 ***E**-field* 下定义的弹性张量，$C_{jk}^D$ 是固定 ***D**-field* 下定义的弹性张量，那么他们通过下式相关联：

$$
C_{jk}^{D}=C_{jk}^{E}+e_{\alpha j}(\epsilon)_{\alpha\beta}^{-1}e_{\beta k}
$$

其中 $e_{\alpha j}$ 是离子弛豫压电张量。

# Calculation Details

## 1 Static dielectric properties

### 1.1 Density functional perturbation theory

密度泛函微扰理论 - [`LEPSILON = .True.`](https://www.vasp.at/wiki/index.php/LEPSILON)

`INCAR` 示例：

```
ISMEAR =  0
SIGMA  =  0.01
EDIFF  = 1.E-8
   
## to get the Born effective charges
## and the macroscopic dielectric tensor
LEPSILON = .TRUE.
    
#LRPA = .TRUE.
#LPEAD = .TRUE.
   
## to get the ionic contribution
## to the macroscopic dielectric tensor
#IBRION = 8
   
## As an alternative to LEPSILON = .TRUE.
## you might try the following:
#LCALCEPS = .TRUE.
   
## and:
#IBRION = 6
#NFREE = 2
```

> - [`LRPA`](https://www.vasp.at/wiki/index.php/LRPA) 标签
>
> 默认情况下，介电张量通过独立粒子近似（*independent-particle approximation*, IPA）计算。
>
> `OUTCAR` 文件中包含以下几行：
>
> ```
> HEAD OF MICROSCOPIC STATIC DIELECTRIC TENSOR (independent particle, excluding Hartree and local field effects)
> ```
>
> ```
> MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)
> ```
>
> 如果添加
>
> ```
> LRPA = .TRUE.
> ```
>
> 第二段将只包含势的Hartree部分响应的局域场效应，即随机相位近似（*random-phase-approximation*, RPA）。`OUTCAR` 显示：
>
> ```
> MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))
> ```
>
> - [`LPEAD`](https://www.vasp.at/wiki/index.php/LPEAD) 标签
>
> 作为微扰表达式的另一种选择，可以在 `INCAR` 中进行以下设置通过有限差分计算 $|\nabla_{\mathbf{k}}\tilde{u}_{n\mathbf{k}}\rangle$：
>
> ```
> LPEAD = .TRUE.
> ```
>
> 然后通过四阶有限差分stencil计算波函数的晶胞周期部分相对于Bloch矢量的导数。
>
> 通过 `LEPSILON = .TRUE.` 计算静态介电性质使用 `LPEAD = .TRUE.` 进行 *k* 点采样，结果的收敛速度会更快。



> **Born有效电荷张量（$Z_{ij}^*$）**
>
> 粗略地讲，Born有效张量提供了一种**当移动一个原子时随原子有效移动的电荷数量**的测度。定义参考[Berry phases and finite electric fields](https://www.vasp.at/wiki/index.php/Berry_phases_and_finite_electric_fields)。对于 `LEPSILON = .True.`，Born有效电荷张量在 `OUTCAR` 结尾给出：
>
> ```
> BORN EFFECTIVE CHARGES (in e, cummulative output)
> ```
>
> ❕ **Mind**: 该条目仅在 `LRPA = .FALSE.`（默认）时显示，因为**RPA中的Born有效电荷往往是无意义的**。

### 1.2 Response to finite electric fields

第二种计算静态介电性质的方法来自[系统对有限电场的自洽响应](https://www.vasp.at/wiki/index.php/Berry_phases_and_finite_electric_fields)。

`INCAR` 示例：

```
ISMEAR =  0
SIGMA  =  0.01
EDIFF  = 1.E-8
    
LCALCEPS = .TRUE.
```

### 1.3 Ionic contributions to the static dielectric properties

为了获得离子对静态介电性质的贡献，需要计算**力常数矩阵（总能量相对于离子位置的Hessian)**和**内部应变张量(总能量相对于应变场和离子位置的二阶导数)**。这些性质可以通过有限差分（`IBRION = 5, 6`）或微扰理论（`IBRION = 7, 8`）得到。

`INCAR` 示例：

```
ISMEAR =  0
SIGMA  =  0.01
EDIFF  = 1.E-8
   
## to get the Born effective charges
## and the macroscopic dielectric tensor
LEPSILON = .TRUE.
LPEAD = .TRUE.
    
## to get the ionic contribution
## to the macroscopic dielectric tensor
IBRION = 8
```

在 `OUTCAR` 文件中搜索

```
MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION
```

```
ELASTIC MODULI IONIC CONTR (kBar)
```

```
PIEZOELECTRIC TENSOR IONIC CONTR  for field in x, y, z        (C/m^2)
```

## 2 Frequency dependent dielectric response

频率依赖的介电函数可以在不同的近似水平上计算：

- 独立粒子近似（IPA）
- 随机相位近似中包含局域场效应
- 密度泛函理论中包含局域场效应

无论后续在介电响应计算方面如何选择，首先都必须从标准DFT（或杂化泛函）计算开始。

`INCAR` - DFT:

```
ISMEAR =  0
SIGMA  =  0.01
EDIFF  = 1.E-8
```

❕ **Mind**: 保留 `WAVECAR` 文件用于后续操作。

### 2.1 The independent-particle picture

为计算独立粒子情况下频率依赖的介电函数，需要从之前DFT计算得到的 `WAVECAR` 重新开始。

`INCAR` - LOPTICS:

```
ALGO = Exact
NBANDS  = 
LOPTICS = .TRUE. ; CSHIFT = 0.100
NEDOS = 2000
   
## and you might try with the following
#LPEAD = .TRUE.
   
ISMEAR =  0
SIGMA  =  0.01
EDIFF  = 1.E-8
```

频率依赖的介电函数被写入 `OUTCAR` 文件：

```
frequency dependent IMAGINARY DIELECTRIC FUNCTION (independent particle, no local field effects)
```

```
 frequency dependent      REAL DIELECTRIC FUNCTION (independent particle, no local field effects)
```

可以使用 `p4vasp` 可视化频率依赖介电函数的实部和虚部：

```
p4v vasprun.xml
```

或运行以下 `bash` 脚本（`plotoptics2`）：

```
awk 'BEGIN{i=1} /imag/,\
                /\/imag/ \
                 {a[i]=$2 ; b[i]=$3 ; i=i+1} \
     END{for (j=12;j<i-3;j++) print a[j],b[j]}' vasprun.xml > imag.dat

awk 'BEGIN{i=1} /real/,\
                /\/real/ \
                 {a[i]=$2 ; b[i]=$3 ; i=i+1} \
     END{for (j=12;j<i-3;j++) print a[j],b[j]}' vasprun.xml > real.dat

cat >plotfile<<!
# set term postscript enhanced eps colour lw 2 "Helvetica" 20
# set output "optics.eps"
plot [0:25] "imag.dat" using (\$1):(\$2) w lp, "real.dat" using (\$1):(\$2) w lp
!

gnuplot -persist plotfile
```

❕ **Mind**: 保留 `WAVECAR` 和 `WAVEDER` 文件以供后续使用。同时备份 `vasprun.xml` 文件副本。

### 2.2 Including local field effects

为了确定包含局域场效应的频率依赖介电函数，需要来自前述IPA计算（`ALGO = Exact` 和 `LOPTICS = .True.` 和**充足的虚拟轨道**）的 `WAVECAR` 和 `WAVEDER` 文件。

`INCAR` - CHI:

```
# Frequency dependent dielectric tensor with and
# without local field effects in RPA
# N.B.: beware one first has to have done a
# calculation with ALGO=Exact, LOPTICS=.TRUE.
# and a reasonable number of virtual states (see above)
ALGO = CHI
       
# be sure to take the same number of bands as for
# the LOPTICS=.TRUE. calculation, otherwise the
# WAVEDER file is not read correctly
NBANDS = 
   
ISMEAR =  0
SIGMA  =  0.01
EDIFF  = 1.E-8
     
LWAVE = .FALSE.
LCHARG= .FALSE.
```

关于独立粒子情况下的介电函数信息显示在 `OUTCAR` 文件中：

```
HEAD OF MICROSCOPIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE)
```

默认情况下，对于 `ALGO = CHI`，包含的局域场效应在RPA水平（`LRPA = .True.`），即仅限于Hartree贡献：

```
INVERSE MACROSCOPIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))
```

为了包含RPA之外的局域场效应，即来自DFT交换和关联的贡献，必须在 `INCAR` 文件中指定 `LRPA = .False.`。此时在 `OUTCAR` 文件输出

```
INVERSE MACROSCOPIC DIELECTRIC TENSOR (test charge-test charge, local field effects in DFT)
```

下述 `bash` 脚本（`plotchi`）使用 `awk` 提取频率依赖的介电常数，既包括独立粒子情况，也包括局域场效应（DFT或RPA），并使用 `gnuplot` 绘制实部分量和虚部分量：

```
awk 'BEGIN{i=1} /HEAD OF MICRO/,\
                /XI_LOCAL/ \
                 {if ($4=="dielectric") {a[i]=$1 ; b[i]=$2 ; c[i]=$3 ; i=i+1}} \
     END{for (j=1;j<i;j++) print a[j],b[j],c[j]}' OUTCAR > chi0.dat

awk 'BEGIN{i=1} /INVERSE MACRO/,\
                /XI_TO_W/ \
                 {if ($4=="dielectric") {a[i]=$1 ; b[i]=$2 ; c[i]=$3 ; i=i+1}} \
     END{for (j=1;j<i;j++) print a[j],b[j],c[j]}' OUTCAR > chi.dat
cat >plotfile<<!
# set term postscript enhanced eps colour lw 2 "Helvetica" 20
# set output "optics.eps"

plot "chi0.dat" using (\$1):(\$2)  w lp lt -1 lw 2 pt 4 title "chi0 real", \
     "chi0.dat" using (\$1):(-\$3) w lp lt  0 lw 2 pt 4 title "chi0 imag", \
     "chi.dat"  using (\$1):(\$2)  w lp lt  1 lw 2 pt 2 title "chi  real", \
     "chi.dat"  using (\$1):(-\$3) w lp lt  0 lw 2 pt 2 lc 1 title "chi  imag"
!

gnuplot -persist plotfile
```

如果之前保存了运行 `LOPTICS = .True.`的 `vasprun.xml` 副本（`vasprun_loptics.xml`），则可以执行 `plotall` 来比较使用 `LOPTICS = .True.` 和 `ALGO = CHI` 计算得到的介电函数：

```
vasprun_LOPTICS=vasprun_loptics.xml
OUTCAR_CHI=OUTCAR

awk 'BEGIN{i=1} /imag/,\
                /\/imag/ \
                 {a[i]=$2 ; b[i]=$3 ; i=i+1} \
     END{for (j=12;j<i-3;j++) print a[j],b[j]}' $vasprun_LOPTICS > imag.dat

awk 'BEGIN{i=1} /real/,\
                /\/real/ \
                 {a[i]=$2 ; b[i]=$3 ; i=i+1} \
     END{for (j=12;j<i-3;j++) print a[j],b[j]}' $vasprun_LOPTICS > real.dat

awk 'BEGIN{i=1} /HEAD OF MICRO/,\
                /XI_LOCAL/ \
                 {if ($4=="dielectric") {a[i]=$1 ; b[i]=$2 ; c[i]=$3 ; i=i+1}} \
     END{for (j=1;j<i;j++) print a[j],b[j],c[j]}' $OUTCAR_CHI > chi0.dat

awk 'BEGIN{i=1} /INVERSE MACRO/,\
                /XI_TO_W/ \
                 {if ($4=="dielectric") {a[i]=$1 ; b[i]=$2 ; c[i]=$3 ; i=i+1}} \
     END{for (j=1;j<i;j++) print a[j],b[j],c[j]}' $OUTCAR_CHI > chi.dat

cat >plotfile<<!
# set term postscript enhanced eps colour lw 2 "Helvetica" 20
# set output "optics.eps"

plot "chi0.dat" using (\$1):(\$2)  w lp lt -1 lw 2 pt 4 title "chi0 real", \
     "chi0.dat" using (\$1):(-\$3) w lp lt  0 lw 2 pt 4 title "chi0 imag", \
     "chi.dat"  using (\$1):(\$2)  w lp lt  1 lw 2 pt 2 title "chi  real", \
     "chi.dat"  using (\$1):(-\$3) w lp lt  0 lw 2 pt 2 lc 1 title "chi  imag", \
     "real.dat"  using (\$1):(\$2) w l lt -1  title "optics  real", \
     "imag.dat"  using (\$1):(-\$2) w l lt  0 lc -1 title "optics  imag"
!

gnuplot -persist plotfile
```

### 2.3 Ionic contributions to the frequency dependent dielectric function

如果要获得完整的（电子 & 离子）频率依赖的介电函数，必须分别进行计算：

1. 执行标准 DFT 计算；

2. 计算***电子***对频率依赖介电函数的贡献；

3. 计算***离子***对频率依赖介电函数的贡献；

4. 将两个介电函数相加：
   $$
   \varepsilon(\omega)=\varepsilon_{\mathrm{elec}}(\omega)+\varepsilon_{\mathrm{ion}}(\omega)
   $$

❕ **Mind**: 该方法可以用于半导体和绝缘体，但不适用于金属。

[`INCAR`](https://www.vasp.at/wiki/index.php/INCAR) - static calculation:

```
PREC = High
ISMEAR = 0 ; SIGMA = 0.01
EDIFF = 1.E-8
GGA = PS
```

[`INCAR`](https://www.vasp.at/wiki/index.php/INCAR) - electronic contributions to the frequency dependent dielectric function:

```
PREC = High
ISMEAR = 0 ; SIGMA = 0.01
EDIFF = 1.E-8
GGA = PS
ALGO = Exact
LOPTICS = .TRUE.
```

[`INCAR`](https://www.vasp.at/wiki/index.php/INCAR) - ionic contributions to the frequency dependent dielectric function:

```
PREC = High
ISMEAR = 0 ; SIGMA = 0.01
EDIFF = 1.E-8
GGA = PS
EDIFF = 1.E-8
#The ionic dielectric function can be calculated in two ways:################
#1# DFPT (faster), but does not allow for METAGGA use. ######################
IBRION = 8; LEPSILON=.TRUE.
#2# Finite differences (slower). ############################################
#IBRION = 6; LPEAD=.TRUE; LCALCEPS=.TRUE.
#NFREE = 2 ; POTIM = 0.015
#In both 1 and 2 the calculated dielectric function is in vasprun.xml #######
```

***Notes***:

介电函数的电子贡献和离子贡献的频率单位不同，分别为 *eV* 和 $2\pi THz$。高频 *ion-clamped* 介电常数（$\epsilon_{\infty}$）通过 `LOPTICS` 计算得到；离子晶格的静态场（$\omega=0$）对介电常数的贡献（$\varepsilon_{\mathrm{ion}}$）通过DFPT或有限差分法计算。两种情况下的声子频率（即动力学矩阵的特征值）和Born有效电荷都会被计算并写入 `OUTCAR` 文件。如果一个声子模式是 *dipole-active* 的，即参与该模式的所有离子的偶极矩大小都会在一个周期内发生变化，那么它将出现在频率依赖的介电函数中。正如动力学矩阵的特征向量一样。总的静态介电常数是两部分贡献之和：
$$
\varepsilon _{{0}}=\varepsilon _{\infty}+\varepsilon _{\rm {ion}}
$$






***Ref***:

- [Dielectric properties - VASP Wiki](https://www.vasp.at/wiki/index.php/Category:Dielectric_properties)
- [Optical properties and dielectric response - Tutorial](https://www.vasp.at/wiki/index.php/Optical_properties_and_dielectric_response_-_Tutorial)
- [*Phys. Rev. B* **73**, 045112 (2006)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.73.045112)
- [*Phys. Rev. B* **72**, 035105 (2005)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.72.035105)
- [*Sci Rep* **6**, 28618 (2016)](https://www.nature.com/articles/srep28618)