---
title: 'Elastic Mechanics'
date: 2024-04-06
tags:
  - VASP
---

> 二阶弹性常数（SOECs）控制着材料的力学性质，对材料的**稳定性**和刚度有重要影响。

# 弹性常数（Elastic Constants）

在材料**线弹性范围**内，固体对外加应变$\varepsilon=(\varepsilon_1,\varepsilon_2,\varepsilon_3,\varepsilon_4,\varepsilon_5,\varepsilon_6)$的应力响应$\sigma=\sigma_1,\sigma_2,\sigma_3,\sigma_4,\sigma_5,\sigma_6)$满足广义[胡克定律(Hooke’s law)](https://en.wikipedia.org/wiki/Hooke%27s_law)：

$$
\sigma_i=\sum_{j=1}^6C_{ij}\varepsilon_j
$$

其中应变$\sigma_i$和应力$\varepsilon_j$分别表示为具有6个独立分量的向量，即$1\leq i, j\leq6$。

**$C_{ij}$是以$GPa$为单位的$6 \times 6$对称矩阵表示的二阶弹性刚度张量**，可由*应力-应变曲线的一阶导数*确定。

$$
C_{ij}=\begin{pmatrix}c_{1111}&c_{1122}&c_{1133}&c_{1123}&c_{1131}&c_{1112}\\c_{2211}&c_{2222}&c_{2223}&c_{2223}&c_{2231}&c_{2212}\\c_{3311}&c_{3322}&c_{3333}&c_{3332}&c_{3331}&c_{3312}\\c_{2311}&c_{2322}&c_{2333}&c_{2323}&c_{2331}&c_{2312}\\c_{3111}&c_{3122}&c_{3133}&c_{3123}&c_{3131}&c_{3112}\\c_{1211}&c_{1222}&c_{1233}&c_{1223}&c_{1231}&c_{1212}\end{pmatrix}\equiv\begin{pmatrix}C_{11}&C_{12}&C_{13}&C_{14}&C_{15}&C_{16}\\C_{12}&C_{22}&C_{23}&C_{24}&C_{25}&C_{26}\\C_{13}&C_{23}&C_{33}&C_{34}&C_{35}&C_{36}\\C_{14}&C_{24}&C_{34}&C_{44}&C_{45}&C_{46}\\C_{15}&C_{25}&C_{35}&C_{45}&C_{55}&C_{56}\\C_{16}&C_{26}&C_{36}&C_{46}&C_{56}&C_{66}\end{pmatrix}
$$

弹性常数$C_{ij}$是应力对应变展开线性项的系数，$i$代表应力方向，$j$代表应力引起的应变方向。

| XX  | YY  | ZZ  | YZ(ZY) | ZX)(XZ) | XY(YX) |
|:---:|:---:|:---:|:------:|:-------:|:------:|
| 1   | 2   | 3   | 4      | 5       | 6      |

可以证明：$C_{ij}=C_{ji}$，因此独立弹性常数最多有**21**个。**独立弹性常数的个数取决于晶体的对称性。对称性越低，独立弹性常数越多。**

$S_{ij}$为**柔度张量**的分量，对应于弹性张量的逆矩阵，即$[ S_{ij} ] = [ C_{ij} ]^{-1}$，单位为$GPa^{-1}$。

# Born弹性稳定性判据

Sufficient and necessary conditions for mechanical stability:

## Cubic

$$
C_{11}-C_{12}>0,\quad C_{11}+2C_{12}>0,\quad C_{44}>0
$$

## Hexagonal & Tetragonal(Ⅰ)

$$
\begin{aligned}&C_{11}>|C_{12}|,\quad2C_{13}^2<C_{33}(C_{11}+C_{12}),\\&C_{44}>0,\quad C_{66}>0\end{aligned}
$$

## Tetragonal(Ⅱ)

$$
C_{11}>|C_{12}|,\quad2C_{13}^2<C_{33}(C_{11}+C_{12}),\\C_{44}>0,\quad2C_{16}^2<C_{66}(C_{11}-C_{12})
$$

## Rhombohedral(Ⅰ)

$$
\begin{aligned}&C_{11}>|C_{12}|,\quad C_{44}>0,\\&C_{13}^{2}<\frac12C_{33}(C_{11}+C_{12}),\\&C_{14}^{2}<\frac12C_{44}(C_{11}-C_{12})\equiv C_{44}C_{66}\end{aligned}
$$

## Rhombohedral(Ⅱ)

$$
\begin{aligned}
C_{11}& >|C_{12}|,\quad C_{44}>0, \\
C_{13}^2& <\frac12C_{33}(C_{11}+C_{12}), \\
C_{14}^2+C_{15}^2& <\frac12C_{44}(C_{11}-C_{12})\equiv C_{44}C_{66}
\end{aligned}
$$

## Orthorhombic

$$
\begin{aligned}&C_{11}>0,\quad C_{11}C_{22}>C_{12}^2,\\&C_{11}C_{22}C_{33}+2C_{12}C_{13}C_{23}-C_{11}C_{23}^2-C_{22}C_{13}^2-C_{33}C_{12}^2>0,\\&C_{44}>0,\quad C_{55}>0,\quad C_{66}>0\end{aligned}
$$

## Monoclinic & Triclinic

Not shown here given the complexity of the equations and solution.

> You can refer to [this article](http://www.gywlxb.cn/cn/article/doi/10.11858/gywlxb.20220575) for details.

# [弹性模量（Elastic modulus）](https://en.wikipedia.org/wiki/Elastic_modulus)

## 体积模量 & 剪切模量

多晶材料的晶粒是随机取向的，在统计意义上，这类材料可以被认为是（准）各向同性。因此，**体弹模量K和剪切模量G一般通过对单晶弹性常数取平均得到**。

> 体积模量：材料的体积与各项均压的比值
> 
> 剪切模量：剪切应力与应变的比值

目前应用最广泛的三种平均方法是：

- **Voigt average**

- **Reuss average**

- **Hill average**

Hill证明了**Voigt和Reuss弹性模量分别是严格的上界和下界**，二者的算术平均值称为**Voigt-Reuss-Hill ( VRH ) average**，它能更好地近似多晶材料的实际弹性行为。

### Voigt average

$$
\left.\left\{\begin{array}{c}9K_V=(C_{11}+C_{22}+C_{33})+2(C_{12}+C_{23}+C_{31})\\15G_V=(C_{11}+C_{22}+C_{33})-(C_{12}+C_{23}+C_{31})+4(C_{44}+C_{55}+C_{66})\end{array}\right.\right.
$$

### Reuss average

$$
\left.\left\{\begin{array}{c}K_R^{-1}=(S_{11}+S_{22}+S_{33})+2(S_{12}+S_{23}+S_{31})\\15G_R^{-1}=4(S_{11}+S_{22}+S_{33})-4(S_{12}+S_{23}+S_{31})+3(S_{44}+S_{55}+S_{66})\end{array}\right.\right.
$$

### Voigt-Reuss-Hill(VRH) average

$$
K=\frac{1}{2}(K_{V}+K_{R})\\~\\G=\frac{1}{2}(G_{V}+G_{R})
$$

## 杨氏模量 & 泊松比

> 杨氏模量：弹性限度内物体应力与应变的比值
> 
> 泊松比：横向正应变与轴向正应变绝对值的比值

$$
\begin{aligned}E&=\frac{9KG}{3K+G}\\~\\\nu&=\frac{3K-2G}{2(3K+G)}\end{aligned}
$$



# VASP计算

VASP可以计算弹性常数，进而得到材料的力学性能。
*关键参数*：

- [IBRION](https://www.vasp.at/wiki/index.php/IBRION)

- [ISIF](https://www.vasp.at/wiki/index.php/ISIF)

- [NFREE](https://www.vasp.at/wiki/index.php/NFREE)



```
IBRON        =    6
ISIF         =    3
NFREE        =    4 or 2
```



INCAR 参考：

```
SYSTEM   =  elastic_constants
ISTART   =  0
ICHARG   =  2

PREC     =  Accurate
ENCUT    =  400
EDIFF    =  1E-6
EDIFFG   = -0.01

IBRION   =  6
ISIF     =  3
NFREE    =  2
POTIM    =  0.015

NSW      =  1
NELM     =  100

ISMEAR   =  0
SIGMA    =  0.05

KSPACING =  0.15

LCHARG   = .FALSE.
LWAVE    = .FALSE.
```



# Appendix

| Crystal system | Point groups                                             | Space groups | Number of independent SOECs $C_{ij}$ |
|:--------------:|:--------------------------------------------------------:|:------------:|:------------------------------------:|
| Triclinic      | 2($1$, $\bar1$)                                          | 2(1-2)       | 21                                   |
| Monoclinic     | 3($2$, $m$, $2/m$)                                       | 13(3-15)     | 13                                   |
| Orthorhombic   | 3($222$, $mm2$, $mmm$)                                   | 59(16-74)    | 9                                    |
| Tetragonal (Ⅱ) | 3($4$, $\bar4$, $4/m$)                                   | 14(75-88)    | 7                                    |
| Tetragonal (Ⅰ) | 4($422$, $4mm$, $\bar42m$, $4/mmm$)                      | 54(89-142)   | 6                                    |
| Trigonal (Ⅱ)   | 2($3$, $\bar3$)                                          | 6(143-148)   | 7                                    |
| Trigonal (Ⅰ)   | 3($32$, $3m$, $\bar3m$)                                  | 19(149-167)  | 6                                    |
| Hexagonal      | 7($6$, $\bar6$, $6/m$, $622$, $6mm$, $\bar6m2$, $6/mmm$) | 27(168-194)  | 5                                    |
| Cubic          | 5($23$, $m\bar3$, $432$, $\bar43m$, $m\bar3m$)           | 36(195-230)  | 3                                    |

## Triclinic: 21 independent elastic constants

$$
C_{ij}=\begin{pmatrix}C_{11}&C_{12}&C_{13}&C_{14}&C_{15}&C_{16}\\C_{12}&C_{22}&C_{23}&C_{24}&C_{25}&C_{26}\\C_{13}&C_{23}&C_{33}&C_{34}&C_{35}&C_{36}\\C_{14}&C_{24}&C_{34}&C_{44}&C_{45}&C_{46}\\C_{15}&C_{25}&C_{35}&C_{45}&C_{55}&C_{56}\\C_{16}&C_{26}&C_{36}&C_{46}&C_{56}&C_{66}\end{pmatrix}
$$

## Monoclinic: 13 independent elastic constants

$$
C_{ij}=\begin{pmatrix}C_{11}&C_{12}&C_{13}&0&C_{15}&0\\C_{12}&C_{22}&C_{23}&0&C_{25}&0\\C_{13}&C_{23}&C_{33}&0&C_{35}&0\\0&0&0&C_{44}&0&C_{46}\\C_{15}&C_{25}&C_{35}&0&C_{55}&0\\0&0&0&0&C_{46}&C_{66}\end{pmatrix}
$$

## Orthorhombic: 9 independent elastic constants

$$
C_{ij}=\begin{pmatrix}C_{11}&C_{12}&C_{13}&0&0&0\\C_{12}&C_{22}&C_{23}&0&0&0\\C_{13}&C_{23}&C_{33}&0&0&0\\0&0&0&C_{44}&0&0\\0&0&0&0&C_{55}&0\\0&0&0&0&0&C_{66}\end{pmatrix}
$$

## Tetragonal (Ⅱ): 7 independent elastic constants

$$
C_{ij}=\begin{pmatrix}C_{11}&C_{12}&C_{13}&0&0&C_{16}\\C_{12}&C_{11}&C_{13}&0&0&-C_{16}\\C_{13}&C_{13}&C_{33}&0&0&0\\0&0&0&C_{44}&0&0\\0&0&0&0&C_{44}&0\\C_{16}&-C_{16}&0&0&0&C_{66}\end{pmatrix}
$$

## Tetragonal (Ⅰ): 6 independent elastic constants

$$
C_{ij}=\begin{pmatrix}C_{11}&C_{12}&C_{13}&0&0&0\\C_{12}&C_{11}&C_{13}&0&0&0\\C_{13}&C_{13}&C_{33}&0&0&0\\0&0&0&C_{44}&0&0\\0&0&0&0&C_{44}&0\\0&0&0&0&0&C_{66}\end{pmatrix}
$$

## Trigonal (Ⅱ): 7 independent elastic constants

$$
C_{ij}=\begin{pmatrix}C_{11}&C_{12}&C_{13}&C_{14}&C_{15}&0\\C_{12}&C_{11}&C_{13}&-C_{14}&-C_{15}&0\\C_{13}&C_{13}&C_{33}&0&0&0\\C_{14}&-C_{14}&0&C_{44}&0&-C_{15}\\C_{15}&-C_{15}&0&0&C_{44}&C_{14}\\0&0&0&-C_{15}&C_{14}&\frac{C_{11}-C_{12}}2\end{pmatrix}
$$

## Trigonal (Ⅰ): 6 independent elastic constants

$$
C_{ij}=\begin{pmatrix}C_{11}&C_{12}&C_{13}&C_{14}&0&0\\C_{12}&C_{11}&C_{13}&-C_{14}&0&0\\C_{13}&C_{13}&C_{33}&0&0&0\\C_{14}&-C_{14}&0&C_{44}&0&0\\0&0&0&0&C_{44}&C_{14}\\0&0&0&0&C_{14}&\frac{C_{11}-C_{12}}2\end{pmatrix}
$$

## Hexagonal: 5 independent elastic constants

$$
C_{ij}=\begin{pmatrix}C_{11}&C_{12}&C_{13}&0&0&0\\C_{12}&C_{11}&C_{13}&0&0&0\\C_{13}&C_{13}&C_{33}&0&0&0\\0&0&0&C_{44}&0&0\\0&0&0&0&C_{44}&0\\0&0&0&0&0&\frac{C_{11}-C_{12}}2\end{pmatrix}
$$

## Cubic: 3 independent elastic constants

$$
C_{ij}=\begin{pmatrix}C_{11}&C_{12}&C_{12}&0&0&0\\C_{12}&C_{11}&C_{12}&0&0&0\\C_{12}&C_{12}&C_{11}&0&0&0\\0&0&0&C_{44}&0&0\\0&0&0&0&C_{44}&0\\0&0&0&0&0&C_{44}\end{pmatrix}
$$

***Ref***:

1. [Physical Properties of Crystals](https://onlinelibrary.wiley.com/doi/book/10.1002/9783527621156)

2. [*Phys. Rev. Lett.* **50**, 697 (1983) - First-Principles Calculation of Stress](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.50.697)

3. [*Phys. Rev. B* **90**, 224104 (2014) - Necessary and sufficient elastic stability conditions in various crystal systems](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.224104)

4. [*Chinese Journal of High Pressure Physics*, 2022, 36(5): 051101](http://www.gywlxb.cn/cn/article/doi/10.11858/gywlxb.20220575)

5. [*Phys. Soc. A* **65** 349 - The Elastic Behaviour of a Crystalline Aggregate](https://iopscience.iop.org/article/10.1088/0370-1298/65/5/307)
