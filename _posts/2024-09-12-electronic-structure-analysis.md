---
title: 'Electronic Structure Analysis (VASP)'
date: 2024-9-12
tags:
  - VASP
  - Electronic Structure
---

# Electronic Structure Analysis (VASP)

- Qualitative analysis: ***Charge density difference***

- Quantitative analysis: ***Population analysis***
  
  $$
  \nabla\rho(\mathbf{r})={\frac{\partial\rho(\mathbf{r})}{\partial\mathbf{x}}}u_{x}+{\frac{\partial\rho(\mathbf{r})}{\partial\mathbf{y}}}u_{y}+{\frac{\partial\rho(\mathbf{r})}{\partial\mathbf{z}}}u_{z}=0
  $$
  
  - Bader
  
  - Lodwin
  
  - Mulliken
  
  - Hirshfeld
  
  - ...

## 1 Charge Density Difference

差分电荷密度分析可以获取体系中原子的成键情况与电荷的空间分布，还可以得到成键电子耦合过程的电荷移动及成键极化方向等性质。

### 1.1 Definition

#### 1.1.1 Difference Charge Density (Set Charge Difference)

$$
\Delta\rho=\rho_{AB}-\rho_{A}-\rho_{B}
$$

即生成物（或体系整体$AB$）与反应物（局部片段$A$和$B$）电荷密度之差。该方法常见于**界面研究**（如催化吸附取代等）相关领域。

> ***Note:*** 
> 
> - 在计算后两个量时，原子的位置应该和它们在AB体系中的位置一样；
> 
> - 计算三个量的FFT-grid (NGXF$\times$NGYF$\times$NGZF)需要保持一致。

#### 1.1.2 Deformation Charge Density (Atomic Charge Difference)

$$
\Delta\rho=\rho(AB)_{self-consistent}-\rho(AB)_{atom}
$$

即成键后体系与对应孤立原子电荷密度之差。

### 1.2 Input

首先对$AB$进行结构优化获取优化后结构`CONTCAR`。

#### 1.2.1 Difference Charge Density

$$
\Delta\rho=\rho_{AB}-\rho_{A}-\rho_{B}
$$

三个结构通过以下`INCAR`参数分别进行**静态自洽**计算，得到对应`CHGCAR`。

```
# key parameters

IBRION   = -1
NSW      =  0

NGXF
NGYF
NGZF

LCHARG   = .TRUE.
```

#### 1.2.2 Deformation Charge Density

$$
\Delta\rho=\rho(AB)_{self-consistent}-\rho(AB)_{atom}
$$

$\rho(AB)_{self-consistent}$通过`1.2.1`中的参数得到，$\rho(AB)_{atom}$通过以下参数得到。

```
# key parameters

ICHARG   =  12
NELM     =  0
```

通过以上参数计算得到的`CHGCAR`是由孤立原子电荷进行简单叠加构成。

### 1.3 Analysis

#### 1.3.1 [VTST Scripts](https://theory.cm.utexas.edu/vtsttools/scripts.html)

> ***`chgsum.pl`***
> 
> - usage: `chgsum.pl (CHGCAR1) (CHGCAR2) (fact1) (fact2)`
> 
> - output: `CHGCAR_sum`
>   
>   `CHGCAR_sum`的值为(`CHGCAR1` * `fact1` + `CHGCAR2` * `fact2`)。默认情况下，`fact1` = `fact2` =1.0，因此输出为输入电荷密度文件的总和：
>   
>   $$
>   \rho_{(CHGCAR\_sum)}=\rho_A+\rho_B
>   $$
> 
> ***`chgdiff.pl`***
> 
> - usage: `chgdiff.pl (CHGCAR1) (CHGCAR2)`
> 
> - output: `CHGCAR_diff`
>   
>   生成CHGCAR之差：
>   
>   $$
>   \rho_{(CHGCAR\_diff)}=\rho_{CHGCAR\color{red}2}-\rho_{CHGCAR\color{red}1}
>   $$

```bash
chgsum.pl A/CHGCAR B/CHGCAR
chgdiff.pl /CHGCAR CHGCAR_sum
```

#### 1.3.2 [VASPKIT](https://vaspkit.com/)

```bash
echo "314\nAB/CHGCAR A/CHGCAR B/CHGCAR\n" | vaspkit
```

生成文件`CHGDIFF.vasp`。

最后通过可视化软件作图分析。

## 2 Bader Charge Analysis

Bader电荷布居可以定量分析体系中**原子的电荷分布与转移情况**。

### 2.1 Definition

Richard Bader 对原子的定义完全基于电子电荷密度，他使用电荷密度的<mark>零通量面（zero flux surfaces）</mark>划分原子之间的界限。零通量面是一个二维曲面，其上的电荷密度在垂直于该曲面为最小值（梯度为0）。在分子系统中，原子间电荷密度达到最小时被认为是原子相互分离的自然位置。将对原子核被零通量面包围的区域积分得到的电子数量统计为该原子所带的电子数量，这种统计方法就是**Bader电荷布局**。

一般认为，电子转移较多，则相互作用越明显，然而反之则不一定成立。

### 2.2 Input

首先进行结构优化，然后进行静态自洽计算。

静态计算`INCAR`参数：

```
# key parameters

IBRION   = -1
NSW      =  0

LAECHG   = .TRUE.
LCHARG   = .TRUE.
```

通过[`LAECHG`](https://www.vasp.at/wiki/index.php/LAECHG)标签显式重构*全电子*电荷密度并写入文件，因为VASP代码中的电荷密度（CHGCAR） 文件仅包含价电荷密度。

对于`LAECHG = .TRUE.`，VASP将在所谓的“精细”FFT网格上重构**三个**不同的*全电子*密度：

1. `AECCAR0`: the core density.

2. `AECCAR1`: the proto-atomic valence density (overlapping atomic charge densities).

3. `AECCAR2`: the self-consistent valence density.

### 2.3 Analysis

```bash
bader chargefile
```

> ⚠️**Note** for VASP users:
> 
> **VASP代码中的电荷密度（CHGCAR）文件的一个主要问题是它们只包含价电荷密度。**Bader分析假设电荷密度最大值位于原子（或伪原子）中心。赝势移除了原子中心的电荷，其计算费用昂贵，而且与原子的重要成键性质无关。
> 
> VASP包含一个模块(`aedens`)，它允许从PAW计算中写出核电荷。通过在`INCAR`文件中添加`LAECHG= .TRUE.`，核电荷被写入`AECCAR0`，价电荷被写入`AECCAR2`。这两个电荷密度文件可以使用`chgsum.pl`脚本进行加和:
> 
> ```bash
> chgsum.pl AECCAR0 AECCAR2
> ```
> 
> 总电荷写入`CHGCAR_sum`。
> 
> 然后对总电荷文件进行bader分析：
> 
> ```bash
> bader CHGCAR -ref CHGCAR_sum
> ```
> 
> 最后需要注意的是，需要一个精细的FFT网格来准确地再现正确的总核电荷。有必要做一些计算，增加NG(X,Y,Z)F，直到总电荷是正确的。

**Output files：**

- `ACF.dat`: 包含每个原子的坐标，与之相关联的电荷（根据Bader划分）， 占总体的比例（根据Bader划分）和到达零通量面的最小距离。如果使用了赝势，则应将此距离与核区域的最大截止半径进行比较。

- `BCF.dat`: 包含每个Bader最大值的坐标，该体积内的电荷，最近的原子以及到该原子的距离。

- `AtomVolumes.dat`: 包含分配给每个原子的体积的数量。数字对应`BvAtxxxx.dat`文件的数量。

电子得失情况可以通过`POTCAR`中`ZVAL`与`ACF.dat`对应元素电荷相减得到。

## 3 Electron Localization Function

电子局域函数体现了以给定点电子为参考，在其邻域内找到与其具有相同自旋的电子的概率，它可以表征参考电子的空间局域化程度，也可以描述多电子体系中的电子对概率。

### 3.1 Definition

$$
\mathrm{ELF}=(1+\chi_{\sigma}^{2})^{-1}
$$

其中，

$$
\chi_{\sigma}=D_{\sigma}/D_{\sigma}^{0}
$$

而

$$
D_{ \sigma}^{ 0}=\frac{3}{5}(6\pi^{2})^{2/3}\rho_{\sigma}^{5/3}
$$

易知

$$
0\leqslant ELF\leqslant1
$$

- 若$ELF=0.5$，则$\chi_\sigma=1$，$D_\sigma=D_{\sigma}^{0}$；说明此处电子行为与均匀电子气行为类似，对应金属体系中的电子行为。

- 若$ELF=1$，则$\chi_\sigma=0$，$D_\sigma=0$；说明电子行为完全局域化。

- 若$ELF=0$，则$\chi_\sigma=\infty$，$D_\sigma=\infty$；说明电子行为完全离域化。

### 3.2 Input

首先进行结构优化，然后进行静态自洽计算。

静态计算`INCAR`参数：

```
# key parameters

PREC   =  Accurate
LELF   = .TRUE.
```

### 3.3 Analysis

**Output file**: `ELFCAR`

最后通过可视化软件作图分析。

***Ref:***

1. [Dr. Renqin Zhang - Charge Density Difference](http://renqinzhang.weebly.com/uploads/9/6/1/9/9619514/charge_density_difference.pdf) Get PDF [here](../files/cms/charge_density_difference.pdf)

2. [*Comput. Mater. Sci.* **36**, 354-360 (2006)](https://www.sciencedirect.com/science/article/pii/S0927025605001849)

3. [*J. Chem. Phys.* **92**, 5397–5403 (1990)](https://pubs.aip.org/aip/jcp/article/92/9/5397/95206/A-simple-measure-of-electron-localization-in)
