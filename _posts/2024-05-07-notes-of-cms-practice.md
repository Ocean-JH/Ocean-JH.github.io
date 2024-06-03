---
title: 'Notes of CMS practice'
date: 2024-05-07
tags:
  - notes
  - Computational Materials Science
  - Molecular Dynamics
  - LAMMPS
---

<div style="color:black; background-color:#FFF3E9; border: 1px solid #FFE0C3; border-radius: 10px; margin-bottom:0rem">
    <p style="margin:1rem; padding-left: 1rem; line-height: 2.5;">
        ©️ <b><i>Copyright 2024 @ Jianghai Wang</i></b><br/>
        <i>Author：
            <b>
            <a href="mailto:your_address@email.com">Jianghai Wang 📨 </a>
            </b>
        </i>
        <br/>
        <i>Date：2024-05-07</i><br/>
        <i>License：<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">CC BY-NC-SA 4.0</a></i><br/>
    </p>
</div>

>**算法原理**：
> - [《计算材料学》（分子动力学）算法原理](https://bohrium.dp.tech/notebooks/52743861357)

# LAMMPS

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/2536258452e64b43ac27974ab3ae2662/49457829-152b-44b2-b115-1c22eb400be6.gif)

> **L**arge-scale **A**tomic/**M**olecular **M**assively **P**arallel **S**imulator
> 
> 
> 
> **Citation**:
> 
> **LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales**, A. P. Thompson, H. M. Aktulga, R. Berger, D. S. Bolintineanu, W. M. Brown, P. S. Crozier, P. J. in 't Veld, A. Kohlmeyer, S. G. Moore, T. D. Nguyen, R. Shan, M. J. Stevens, J. Tranchida, C. Trott, S. J. Plimpton, [*Comp Phys Comm*, **271** (2022) 10817](https://doi.org/10.1016/j.cpc.2021.108171).

- [LAMMPS Website](https://www.lammps.org/)

- [LAMMPS Github](https://github.com/lammps/lammps)

- [LAMMPS Doc](https://docs.lammps.org/)



## [Unit of measurement](https://en.wikipedia.org/wiki/Unit_of_measurement)

### [SI unit](https://en.wikipedia.org/wiki/International_System_of_Units)

| Quantity name             | Unit name | Unit symbol |
|:-------------------------:|:---------:|:-----------:|
| time                      | second    | s           |
| length                    | metre     | m           |
| mass                      | kilogram  | kg          |
| electric current          | ampere    | A           |
| thermodynamic temperature | Kelvin    | K           |
| amount of substance       | mole      | mol         |
| luminous intensity        | candela   | cd          |



### L-J unit

#### 1. 原始Lenard-Jones势

$$
V(r)=4\varepsilon\left[(\frac\sigma r)^{12}-(\frac\sigma r)^6\right]
$$

括号中前一项为Pauli exclusion principle引发的短程排斥作用；后一项为London dispersion导致的远程吸引作用。不同元素具有不同的L-J势参数。



> - [Interatomic Potentials - LAMMPS Tube](https://lammpstube.com/mdpotentials/)
> 
> 不同元素之间的L-J势参数可以由以下方式计算得到*：
> 
> $$
\sigma_{AB}=\frac{(\sigma_{AA}+\sigma_{BB})}{2}\\
~\\
\varepsilon_{AB}=\sqrt{\varepsilon_{AA}\cdot\varepsilon_{BB}}
$$
> 
> **Ref*：ARKUNDATO, ARTOTO; SU'UD, ZAKI; ABDULLAH, MIKRAJUDDIN; and SUTRISNO, WIDAYANI (2013) "Molecular dynamic simulation on iron corrosion-reduction in high temperature molten lead-bismuth eutectic," [*Turkish Journal of Physics*: Vol. 37: No. 1, Article 14](https://journals.tubitak.gov.tr/physics/vol37/iss1/14/).



#### 2. 约化Lenard-Jones势

$$
V^{\prime}(r)=4[\frac1{r^{\prime12}}-\frac1{r^{\prime6}}]\\
~\\
V^{\prime}(r)=V(r)/\varepsilon\quad r_i^{\prime}=r_i/\sigma
$$

其他单位约化标度：

- 质量单位：$m_i^{\prime}=m_i/m$

- 长度单位：$r_i^{\prime}=r_i/\sigma$

- 能量单位：$V_i^{\prime}=V_i/\varepsilon$

- 时间单位：$t_i^{\prime}=t_i/\tau$

### Other Units

| Quantity          | Metal Units        | Real Units            |
|:-----------------:|:------------------:|:---------------------:|
| Length            | Angstroms ($\AA$)  | Angstroms ($\AA$)     |
| Time              | Picoseconds ($ps$) | Femtoseconds ($fs$)   |
| Energy            | $eV$               | $kcal/mole$           |
| Temperature       | Kelvin ($K$)       | Kelvin ($K$)          |
| Pressure          | Bars               | Atm                   |
| Velocity          | $\AA/ps$           | $\AA/fs$              |
| Force             | $eV/\AA$           | $kcal/(mole\cdot\AA)$ |
| Torque            | $eV$               | $kcal/mole$           |
| Dynamic Viscosity | $eV\cdot ps/\AA^3$ | $kcal/(fs\cdot\AA^3)$ |

## Interatomic Potentials
- [*Interatomic Potentials Repository* (nist.gov)](https://www.ctcms.nist.gov/potentials/)

---

# LAMMPS实例

## 1 气体分子扩散

### 1.1 单组分气体


```python
%%writefile in.single_comp_diffusion

#---------------Initialize Simulation -----------------------#
units lj                # 单位
dimension 2             # 维度
boundary p p p          # 边界条件
atom_style atomic       # 原子类型

#-------------- Create Atoms Initial Conditions------------- #
lattice hex 1.0                        # 晶格类型
region box block 0 20 0 10 -0.1 0.1    # 空间区域
create_box 1 box                       # 模拟盒子
region 2 block 5 15 0 10 -0.1 0.1
create_atoms 1 region 2                # 基于晶格点阵创建原子
mass 1 1.0                             # 质量
velocity all create 2.5 87287          # 速度

#---------------- Define Interatomic Potential --------------#
pair_style lj/cut 2.5                  # 原子相互作用势
pair_coeff 1 1 1.0 1.0 2.5             # 对势参数
neighbor 0.3 bin                       # 近邻列表
neigh_modify every 20 delay 0 check no # 近邻算法
fix 1 all nvt temp 0.5 0.5 0.01        # 时间步操作
fix 2 all enforce2d

#--------------- Run MD Simulation --------------------------#
dump 1 all custom 100 toEquil.lammpstrj id type x y z vx vy vz    # 输出分子动力学模拟信息
thermo 500                                                        # 输出热力学信息间隔
run 10000                                                         # 运行指定步数分子动力学模拟
```


```python
# !lmp -i in.single_comp_diffusion
```

**结果可视化**：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/cb2c4534386e4a5f8b303a1e07f7fbd0/b774c372-18b8-4542-98b1-55bca32072ac.gif)

### 1.2 多组分气体


```python
%%writefile in.multi_comp_diffusion

#---------------Initialize Simulation -----------------------#
units lj
dimension 2
boundary p p p
atom_style atomic
variable t equal 0.5             # 定义系统温度

#-------------- Create Atoms Initial Conditions------------- #
lattice sq 1.0
region box block 0 100 0 100 -0.5 0.5
create_box 2 box
create_atoms 1 random 2500 12345 box
create_atoms 2 random 2500 54321 box        # 随机创建原子
mass 1 1.0
mass 2 1.0

#---------------- Define Interatomic Potential --------------#
pair_style hybrid lj/cut 2.5 soft 5.0
pair_coeff 1 1 lj/cut 1.0 1.0 2.5
pair_coeff 2 2 lj/cut 1.0 1.0 2.5
pair_coeff 1 2 soft 5.0

#--------------- Run MD Simulation --------------------------#
compute eng all pe/atom                                        # 计算每个原子势能
compute eatoms all reduce sum c_eng                            # 计算所有原子势能
thermo_style custom step temp epair etotal press c_eatoms      # 定义输出的热力学量
thermo 1000
dump id all atom 100 dump.lammpstrj
minimize 1e-4 1e-6 1000 10000                                  # 防止原子距离过近
velocity all create $t 87287
fix nvt all nvt temp $t $t 0.01
run 50000
```


```python
%%capture
# !lmp -i in.multi_comp_diffusion
```

*在进行动力学演化前，先进行能量最小化，防止随机产生原子
距离过近导致的系统不稳定性。*

**结果可视化**：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/cb2c4534386e4a5f8b303a1e07f7fbd0/a4e5fe9c-3fb6-4ca6-8ee2-063df5628fe0.gif)

## 2 力学性质

### 2.1 平衡晶格常数

晶格是晶体结构的数学表示，晶格中的每个格点代表一个基元。



[**Crystal system**](https://en.wikipedia.org/wiki/Crystal_system)

| Crystal system | Required symmetries of the point group                                       | Point groups | Space groups |
|:--------------:|:----------------------------------------------------------------------------:|:------------:|:------------:|
| Triclinic      | None                                                                         | 2            | 2            |
| Monoclinic     | 1 twofold axis of rotation or 1 mirror plane                                 | 3            | 13           |
| Orthorhombic   | 3 twofold axes of rotation or 1 twofold axis of rotation and 2 mirror planes | 3            | 59           |
| Tetragonal     | 1 fourfold axis of rotation                                                  | 7            | 68           |
| Trigonal       | 1 threefold axis of rotation                                                 | 5            | 25           |
| Hexagonal      | 1 sixfold axis of rotation                                                   | 7            | 27           |
| Cubic          | 4 threefold axes of rotation                                                 | 5            | 36           |



平衡晶格常数对应的晶格具有最小的结合能；结合能曲线最低点对应最小结合能$E(a_0)$和平衡晶格常数$a_0$。



#### 2.1.1 FCC-Ar平衡晶格常数

##### 2.1.1.1 单点计算


```python
%%writefile in.Ar_single_point

units lj
boundary p p p
atom_style atomic
lattice fcc 1.00
region box block 0 1 0 1 0 1
create_box 1 box
create_atoms 1 box
mass 1 1.0
pair_style lj/cut 8.0
pair_coeff 1 1 1.0 1.0 8.0
variable P equal pe
variable L equal (count(all)/1.00)^(1/3)
run 0
print “Lattice constant: $L"
print "Cohesive Energy of Ar: $P"
```

其中Ar的L-J势参数如下：

|     | $\sigma(nm)$ | $\epsilon(J)$ |
|:---:|:------------:|:-------------:|
| Ar  | 1.3405       | 1.6540 E-21    |

*Ref*: [*Dokl Phys Chem* **472**, 16–18 (2017)](https://doi.org/10.1134/S0012501617010043)

##### 2.1.1.2 循环控制




```python
%%writefile in.Ar_loop

units lj
boundary p p p
atom_style atomic

label loop_i
variable i  loop  50
variable x  equal  1.02+0.002*$i

lattice fcc $x
region     box block 0 1 0 1 0 1
create_box 1 box
create_atoms 1 box
mass 1 1.0

pair_style  lj/cut  4.0
pair_coeff  1  1  1.0  1.0  4.0

variable P equal pe
variable L equal (count(all)/$x)^(1/3)

run 0
print  "Cohesive Energy of Ar a = $L E = $P"
clear
next i
jump SELF loop_i
```

**数据提取**：


```python
# Read the file and search for lines containing the pattern "Cohesive Energy of Ar a = * E = *
import re

# Initialize list to hold the matching lines
matching_lines = []

# Define the regex pattern for capturing a and E values
pattern = re.compile(r"Cohesive Energy of Ar a = ([\d\.-]+)\s+E = ([\d\.-]+)")

# Read the file and search for the pattern
with open('./log.lammps', 'r') as file:
    for line in file:
        match = pattern.search(line)
        if match:
            a_value = float(match.group(1))
            E_value = float(match.group(2))
            matching_lines.append((a_value, E_value))

# Write the extracted a and E values to a CSV file
with open('Ar_fcc.csv', 'w') as f:
    f.write("a,E\n")  # Write header
    for a, E in matching_lines:
        f.write(f"{a},{E}\n")
```

**结果可视化**：


```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


filename = "Ar_fcc.csv"
DataFrame = pd.read_csv(filename)
lat = np.array(DataFrame['a'])
ene = np.array(DataFrame['E'])

a, b, c = np.polyfit(lat, ene, 2)
mesh = np.linspace(np.min(lat), np.max(lat), 1000)
opt_lat = -b / (2 * a)
Emin = np.polyval([a, b, c], opt_lat)
print("a,b,c =", a, b, c)
print("Emin =", Emin, "opt_lat =", opt_lat)

plt.scatter(lat, ene, 10)
plt.plot(mesh, a * mesh ** 2 + b * mesh + c)
plt.xlabel('Lattice parameter')
plt.ylabel('Energy')
plt.show()
```

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/cc022afdb749478486df31dfedadb978/bd1a856d-d9e2-4617-9a09-dc978643ad43.png)



- 理论平衡晶格常数：

$$
\frac{\partial u}{\partial r}=-2\varepsilon[nA_{n}(\frac{\sigma^{n}}{r^{n+1}})-mA_{m}(\frac{\sigma^{m}}{r^{m+1}})]\\
~\\
-2\varepsilon[nA_{n}(\frac{\sigma^{n}}{r^{n+1}})-mA_{m}(\frac{\sigma^{m}}{r^{m+1}})]=0\\
~\\
r_{0}=(\frac{2A_{12}}{A_{6}})^{1/6}\sigma=1.09\sigma\\
~\\
\text{a}=\sqrt{2}r_0=1.54\sigma
$$

#### 2.1.2 FCC-Al平衡晶格常数

##### 2.1.2.1 循环控制


```python
%%writefile in.Al_loop

units metal
boundary p p p
atom_style atomic
variable i loop 30
variable x equal 3.90+0.01*$i

# Define a lattice for use by other commands.
lattice fcc $x
region box block 0 1 0 1 0 1
create_box 1 box
create_atoms 1 box

pair_style eam/fs
pair_coeff * * Al_mm.eam.fs Al

variable n equal count(all)
variable P equal pe/$n

run 0
print "Cohesive Energy of Al a = $x E = $P"
clear 
next i
jump SELF
```

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/8f49555cc5cc4477ba30f8a8a44d6639/fe7abbc3-81f8-4cf1-97d7-ea9956a0ef43.png)

##### 2.1.2.2 晶格弛豫

采用共轭梯度法优化晶体结构：


```python
%%writefile in.Al_relax

units metal
dimension 3
boundary p p p
atom_style atomic

lattice fcc 4.00
region box block 0 1 0 1 0 1 units lattice
create_box 1 box
create_atoms 1 box

pair_style eam/fs
pair_coeff * * Al_mm.eam.fs Al

neighbor 2.0 bin
neigh_modify delay 10 check yes

compute eng all pe/atom
compute eatoms all reduce sum c_eng
dump 1 all atom 1 relax.lammpstrj

#--------------Run minimization--------------#
reset_timestep 0
fix 1 all box/relax iso 0.0 vmax 0.001
thermo 10
thermo_style custom step pe lx ly lz press pxx pyy pzz c_eatoms

min_style cg
minimize 1e-25 1e-25 5000 10000

variable natoms equal count(all)
variable teng equal "c_eatoms"
variable a equal lx
variable ecoh equal "v_teng/v_natoms“

print "Total energy (eV)= ${teng};" 
print "Number of atoms = ${natoms};"
print "Lattice constant (Angstroms) = ${a};" 
print "Cohesive energy (eV/atom) = ${ecoh};" 
```

#### 2.1.3 注意事项

- 晶格常数拟合区间过小会放大**截断误差**，影响计算精度；

- 晶格常数拟合区间过大会产生**非谐效应**，使实际曲线偏离二次型。



通过列举常见的晶格类型，拟合得到平衡晶格常数对应的最低能量，即可确定元素的最稳定构型。

### 2.2 体积模量

[体积模量（Bulk modulus）](https://en.wikipedia.org/wiki/Bulk_modulus)是衡量物质**可压缩性**的指标。

#### 2.2.1 计算思路

$T=0$时，压强可以表示为

$$
P=-\frac{\partial U}{\partial V}
$$

其中$U$为体系的总能。定义单位原子的能量和体积：

$$
u=\frac UN\quad v=\frac VN
$$

则

$$
\begin{aligned}B&=-V\left(\frac{\partial P}{\partial V}\right)_T\\&=V\frac{\partial^2U}{\partial V^2}\\&=v\frac{\partial^2u}{\partial v^2}\end{aligned}
$$

#### 2.2.2 解析表达

$$
\begin{aligned}\text{u}&=2\varepsilon[\sum_{i\neq j}(\frac{1}{M_{ij}})^{12}(\frac{\sigma}{r})^{12}-\sum_{i\neq j}(\frac{1}{M_{ij}})^{6}(\frac{\sigma}{r})^{6}]\\&=2\varepsilon[A_{12}\left(\frac\sigma r\right)^{12}-A_6\left(\frac\sigma r\right)^6]\end{aligned}
$$

其中，

$$
A_{12}=\sum_{j}\frac{1}{M_{ij}^{12}}=12.13\\A_{6}=\sum_{j}\frac{1}{M_{ij}^{6}}=14.45
$$

因此，

$$
\begin{aligned}B&=v\frac{\partial^2u}{\partial v^2}\\&=\frac v{\left(\frac{\partial v}{\partial r}\right)^2}\cdot\frac{\partial^2u}{\partial r^2}\\&=\frac{\sqrt{2}}{9r}\cdot\frac{\partial^2u}{\partial r^2}|_{r=r_0}\end{aligned}
$$

> 对于**FCC晶体**，有
> 
> $$
\begin{cases} v=\frac{a^3}4=\frac{r^3}{\sqrt{2}}\\\\ \left(\frac{\partial v}{\partial r}\right)^2=\frac{9r^4}2\end{cases}
$$

$$
\begin{aligned}
u&=2\varepsilon[A_{12}\left(\frac{\sigma}{r}\right)^{12}-A_6\left(\frac{\sigma}{r}\right)^6] \\\\
\frac{\partial u}{\partial r}&=2\varepsilon[\frac{(-12)A_{12}\sigma^{12}}{r^{13}}-\frac{(-6)A_{6}\sigma^{6}}{r^{7}}] \\\\
\frac{\partial^2u}{\partial r^2}&=2\varepsilon[\frac{(-13)(-12)A_{12}\sigma^{12}}{r^{14}}-\frac{(-7)(-6)A_6\sigma^6}{r^8}] \\
\end{aligned}
$$

将

$$
\begin{cases} r_0=1.09\sigma\\\\ \frac{\partial^2u}{\partial r^2}|_{r_0}=523\frac\varepsilon{\sigma^2}\end{cases}
$$

代入

$$
B=\frac{\sqrt{2}}{9r_0}\cdot\frac{\partial^2u}{\partial r^2}|_{r=r_0}
$$

得到理论体积模量

$$
B_0=75\frac\varepsilon{\sigma^3}
$$

#### 2.2.3 L-J势体积模量

- **公式**：$B=v\frac{\partial^2u}{\partial v^2}$

- **方法**：不断调整晶格常数大小，用LAMMPS计算平衡距离附近时，Cu晶体能量与其体积的关系，据此**拟合**体弹性模量数据。


```python
%%writefile in.fcc_bulk_modulus

units lj
boundary p p p
atom_style atomic
label LOOP

# The volume change is controlled by lattice constant
variable i loop 40
variable x equal 1.05+0.002*$i

lattice fcc $x
region box block 0 1 0 1 0 1
create_box 1 box
create_atoms 1 box

pair_style lj/cut 8.50
pair_coeff * * 1.0 1.0
mass 1 1.0

thermo_style custom etotal
variable P equal pe                  # Potential energy per atom
variable v equal (1.0/$x)            # vol in LJ unit is N/x, vol per atom in LJ units is 1/x

dump 1 all custom 10 trj.lammpstrj id type x y z vx vy vz
run 0
print "FCC lattice rho = $v E = $P"
clear
next i
jump SELF LOOP
```

数据处理：


```python
# Read the file and search for lines containing the pattern "FCC lattice rho = * E = *"
import re

# Initialize list to hold the matching lines
matching_lines = []

# Define the regex pattern for capturing a and E values
pattern = re.compile(r"FCC lattice rho = ([\d\.-]+)\s+E = ([\d\.-]+)")

# Read the file and search for the pattern
with open('./log.lammps', 'r') as file:
    for line in file:
        match = pattern.search(line)
        if match:
            v_value = float(match.group(1))
            E_value = float(match.group(2))
            matching_lines.append((v_value, E_value))

# Write the extracted a and E values to a CSV file
with open('fcc_bulk_modulus.csv', 'w') as f:
    f.write("V,E\n")  # Write header
    for a, E in matching_lines:
        f.write(f"{V},{E}\n")
```

##### 2.2.3.1 二次函数拟合


```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


filename = "fcc_bulk_modulus.csv"
DataFrame = pd.read_csv(filename)
vol = np.array(DataFrame['V'])
ene = np.array(DataFrame['E'])

a, b, c = np.polyfit(vol, ene, 2)
mesh = np.linspace(np.min(vol), np.max(vol), 1000)
opt_vol = -b / (2 * a)
Emin = np.polyval([a, b, c], opt_vol)
print("a,b,c =", a, b, c)
print("Emin =", Emin, "opt_vol =", opt_vol)

plt.scatter(vol, ene, 10)
plt.plot(mesh, a * mesh ** 2 + b * mesh + c)
plt.xlabel('Volume')
plt.ylabel('Energy')

plt.show()
```

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/cc022afdb749478486df31dfedadb978/15e2ee08-e5db-45fe-90c0-edfe9c8eff4e.png)

> E-V数据存在一定的**非谐效应**。

##### 2.2.3.2 [Birch-Murnaghan方程](https://en.wikipedia.org/wiki/Birch–Murnaghan_equation_of_state)拟合

$$
E\left(V\right)=E_0+\frac{B_0V}{B_0^{\prime}}{\left(\frac{(V_0/V)^{B_0^{\prime}}}{B_0^{\prime}-1}+1\right)-\frac{B_0V_0}{B_0^{\prime}-1}}
$$

LAMMPS计算结果为$E$、$V$，待拟合参数为$E_0$、$V_0$、$B_0$、$B_0^{\prime}$。


```python
import pandas as pd
import numpy as np


df = pd.read_csv('fcc_bulk_modulus.csv)
v = np.array(df['V'])
e = np.array(df['E'])

a, b, c = np.polyfit(v, e, 2)
v0 = -b / (2 * a)
e0 = a * v0 ** 2 + b * v0 + c
b0 = 2 * a * v0
bP = 3.5
x0 = [e0, b0, bP, v0]


# Birch–Murnaghan equation of state
def Murnaghan(parameters, vol):
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    E = E0 + B0 * vol / BP * (((V0 / vol) ** BP) / (BP - 1) + 1) - V0 * B0 / (BP - 1.)

    return E


# error
def residual(pars, y, x):
    err = y - Murnaghan(pars, x)

    return err


from scipy.optimize import leastsq

murnpars, ier = leastsq(residual, x0, args=(e, v))

print('Bulk Modulus:' + str(murnpars[1]))
print('lattice constant:', murnpars[3] ** (1 / 3))

from matplotlib import pyplot as plt

v_mesh = np.linspace(np.min(v), np.max(v), 1000)
plt.scatter(v, e, 10)
plt.plot(v_mesh, Murnaghan(murnpars, v_mesh))
plt.xlabel('Volume')
plt.ylabel('Energy')

plt.show()
```

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/cc022afdb749478486df31dfedadb978/d02d4341-e97a-4753-a8c4-8106134553df.png)

> **Birch-Murnaghan方程更接近于物理实际**，得到的结果与解析解最为接近，因此**拟合体积模量应优先选择Birch-Murnaghan方程**。

#### 2.2.4 Cu体积模量

Cu的L-J势参数如下：


|     | $\sigma(nm)$ | $\epsilon(J)$ |
|:---:|:------------:|:-------------:|
| Cu  | 0.2338       | 65.5815 E-21  |


*Ref*: [*Nanoscale Res Lett* **6**, 200 (2011)](https://doi.org/10.1186/1556-276X-6-200)



##### 2.2.4.1 L-J势


```python
%%writefile in.Cu_lj_mudulus

units metal
boundary p p p
atom_style atomic
#------------- setup loop -------------#
variable i loop 40
variable x equal 3.40+0.01*$i
lattice fcc $x
region box block 0 1 0 1 0 1
create_box 1 box
create_atoms 1 box
mass 1 64

#------------- Cu LJ parameter -------------#
pair_style lj/cut 10.0
pair_coeff 1 1 0.40933 2.338
variable v equal ($x)^3
variable n equal count(all)
variable P equal pe

#------------- run -------------#
run 0
print "Cohesive Energy of Cu v = $v x = $x E = $P "
clear 
next i
jump SELF
```

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/8f49555cc5cc4477ba30f8a8a44d6639/d643285f-d050-4a30-b4d8-d815a6781e76.png)

> **L-J势函数不能准确描述金属原子的相互作用！**

##### 2.2.4.2 EAM势


```python
%%writefile in.Cu_eam_mudulus

units metal
boundary p p p
atom_style atomic
#------------- setup loop -------------#
variable i loop 40
variable x equal 3.40+0.01*$i
lattice fcc $x
region box block 0 1 0 1 0 1
create_box 1 box
create_atoms 1 box
mass 1 64

#------------- Cu EAM parameter -------------#
pair_style eam
pair_coeff * * Cu_u3.eam
variable v equal ($x)^3
variable n equal count(all)
variable P equal pe

#------------- run -------------#
run 0
print "Cohesive Energy of Cu v = $v x = $x E = $P "
clear 
next i
jump SELF
```

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/8f49555cc5cc4477ba30f8a8a44d6639/a11730cb-202b-47b5-a734-694d47d31d5b.png)

### 2.3 拉伸性能

#### 2.3.1 应力—应变曲线

[**应力**](https://en.wikipedia.org/wiki/Stress_(mechanics))和[**应变**](https://en.wikipedia.org/wiki/Strain_(mechanics))之间的关系称为材料的**应力—应变曲线**。

$$
\begin{aligned}\sigma&=\frac F{A_0}\\\\\epsilon&=\frac{L-L_0}{L_0}=\frac{\Delta L}{L_0}\end{aligned}
$$

**石墨烯拉伸试验**：
- 石墨烯建模


```python
%%writefile in.graphene

# Parameters
variable        Lx      equal   20
variable        Ly      equal   10
# variable      Lz      equal   0.8

# Unit Size
variable        x         equal    ${Lx}
variable        y         equal    ${Ly}
# variable        z         equal    ${Lz}
variable        xbox         equal    round(v_x)
variable        ybox         equal    round(v_y)
# variable        zbox         equal    round(v_z)

# Initialization
units                real
dimension            3
boundary             p p p
neighbor             2.0 bin
neigh_modify         every 10 delay 0 check yes
timestep             0.001
atom_style           charge

# Modeling
lattice custom  1.421 a1  3  0  0   a2   0 1.732  0 a3 0 0 2.357   &
        basis   0        0  0  &
        basis   0.333    0  0  &
        basis   0.5     0.5 0  &
        basis   0.833   0.5 0
region                box block 0 ${xbox} 0 ${ybox} -5.0 5.0
create_box            1 box
region                graphene block 0 ${xbox} 0 ${ybox} -0.1 0.1
create_atoms          1 region graphene
mass                  * 12.011150
write_data            graphene.dat
```

- 石墨烯拉伸


```python
%%writefile in.tensile

# UNITS
units real
timestep 1
variable fpress equal 0.000101325        # atm -> GPa
variable fenergy equal 0.043             # kcal/mole -> eV

# GRAPHENE SHEET
dimension 3
boundary p p f
atom_style charge
read_data graphene.dat
pair_style reax/c NULL checkqeq no
pair_coeff * * ffield.reax.cho C

# TEMPERATURE SETTINGS
variable temp equal 300
variable seed equal 1717
velocity all create ${temp} ${seed} dist uniform

# EQUILIBRATION
fix fnpt all npt temp ${temp} ${temp} 10 x 0 0 500 y 0 0 500
thermo 100
run 1000

# OUTPUT
# assume that the thickness of monolayer graphene is 3.35A.
variable tmp equal lx
variable lx0 equal ${tmp}
variable tmp equal ly
variable ly0 equal ${tmp}
variable Eavg equal etotal/atoms*${fenergy}                # eV/atom
variable pe equal pe/atoms*${fenergy}                      # eV/atom
variable ke equal ke/atoms*${fenergy}                      # eV/atom
variable strainx equal (lx-${lx0})/${lx0}
variable strainy equal (ly-${ly0})/${ly0}
variable stressx equal -pxx*(lz/3.35)*${fpress}            # GPa
variable stressy equal -pyy*(lz/3.35)*${fpress}            # GPa
thermo_style custom step time temp etotal press v_Eavg v_pe v_ke v_strainx v_stressx v_strainy v_stressy

# DEFORMATION
fix boxdeform all deform 1 x scale 2 remap x
fix fnpt all npt temp ${temp} ${temp} 10 y 0 0 500
fix output all ave/time 1 100 100 v_strainx v_stressx v_strainy v_stressy file stress_strain.txt
dump 1 all atom 10 dump.lammpstrj
fix stop all halt 100 v_strainx > 0.7  error  continue
thermo 100
run 10000
```

*石墨烯拉伸演示*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/9eea4c3ad7fd4a9e9a8e2ca50b8c4fad/7cb4df51-0e93-4b89-96ce-22d356209071.gif)

*应力—应变曲线*：
![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/58f11699c79f4f4db5b19d22961bed65/a3d71a42-21a4-4fe5-bf37-7d582f7f0f2a.png)

#### 2.3.2 弹性常数

[**弹性常数**](https://en.wikipedia.org/wiki/Elastic_modulus)表征材料弹性的量，弹性变形时应力与应变满足[*胡克定律*](https://en.wikipedia.org/wiki/Hooke%27s_law)。联系各向异性介质中应力和应变关系的广义弹性张量共有**81**个分量，其中有**21个独立分量**，可以根据**晶体对称性**简化。

$$
\sigma_{ij}=C_{ijkl}\varepsilon_{kl}\\
~\\
c_{ijmn}=c_{jimn}=c_{ijnm}=c_{jinm}=c_{mnij}=c_{nmij}=c_{mnji}=c_{nmji}\\
~\\
\left.\left[\begin{array}{ccccccc}c_{1111}&c_{1122}&c_{1133}&c_{1123}&c_{1113}&c_{1112}\\&c_{2222}&c_{2233}&c_{2223}&c_{2213}&c_{2212}\\&&c_{3333}&c_{3323}&c_{3313}&c_{3312}\\&&&c_{2323}&c_{2313}&c_{2312}\\&&&&c_{1313}&c_{1312}\\&&&&&c_{1212}\end{array}\right.\right]
$$



以FCC-Cu为例，根据对称性与独立弹性常数个数关系，具有立方对称性的材料有3个独立的弹性常数：

$$
\left.C=\left(\begin{array}{cccccc}C_{11}&C_{12}&C_{12}&0&0&0\\C_{12}&C_{11}&C_{12}&0&0&0\\C_{12}&C_{12}&C_{11}&0&0&0\\0&0&0&C_{44}&0&0\\0&0&0&0&C_{44}&0\\0&0&0&0&0&C_{44}\end{array}\right.\right)
$$

体系的弹性内能：

$$
E^\text{elas}/V=\frac12\sum_{ij}C_{ij}\varepsilon_i\varepsilon_j
$$

因此可以构造相应的弹性形变来计算弹性常数。



##### 2.3.2.1 刚度张量

$$
C=\left(\begin{array}{cccccc}C_{11}&C_{12}&C_{12}&0&0&0\\C_{12}&C_{11}&C_{12}&0&0&0\\C_{12}&C_{12}&C_{11}&0&0&0\\0&0&0&C_{44}&0&0\\0&0&0&0&C_{44}&0\\0&0&0&0&0&C_{44}\end{array}\right)
$$



##### 2.3.2.2 柔度张量

$$
[S]=[C]^{-1}
$$



$$
\left.S=\left(\begin{array}{cccccc}\frac{C_{11}+C_{11}}{C_{11}^2+C_{11}C_{12}-2C_{12}^2}&\frac{-C_{12}}{C_{11}^2+C_{11}C_{12}-2C_{12}^2}&\frac{-C_{12}}{C_{11}^2+C_{11}C_{12}-2C_{12}^2}&0&0&0\\\frac{-C_{12}}{C_{11}^2+C_{11}C_{12}-2C_{12}^2}&\frac{C_{11}+C_{11}C_{12}-2C_{12}^2}{C_{11}^2+C_{11}C_{12}-2C_{12}^2}&\frac{-C_{12}}{C_{11}^2+C_{11}C_{12}-2C_{12}^2}&0&0&0\\\frac{-C_{12}}{C_{11}^2+C_{11}C_{12}-2C_{12}^2}&\frac{-C_{11}}{C_{11}^2+C_{11}C_{12}-2C_{12}^2}&\frac{C_{11}+C_{12}}{C_{11}^2+C_{11}C_{12}-2C_{12}^2}&0&0&0\\0&0&0&1/C_{44}&0&0\\0&0&0&0&0&1/C_{44}\end{array}\right.\right)
$$





#### 2.3.3 杨氏模量

弹性材料承受正向应力时会产生正向应变，定义**正向应力与正向应变的比值**为材料的[**杨氏模量**](https://en.wikipedia.org/wiki/Young%27s_modulus)：

$$
Y=\frac\sigma\varepsilon=\frac{F/_{A_0}}{\Delta L/_{L_0}}=\frac{FL_0}{A_0\Delta L}
$$

弹性材料承受单轴拉伸时，唯一的应力为$\sigma_{xx}$

$$
\varepsilon_{xx}=S_{11}\sigma_{xx}\\~\\
Y=\frac{\sigma_{xx}}{\varepsilon_{xx}}=\frac1{S_{11}}=\frac{C_{11}{}^2+C_{11}C_{12}-2C_{12}{}^2}{C_{11}+C_{12}}
$$

其中$S_{11}$为柔性张量的$(1,1)$分量。



#### 2.3.4 泊松比

材料受拉伸或压缩时，其**横向形变量与纵向形变量**的比值称为[**泊松比**](https://en.wikipedia.org/wiki/Poisson%27s_ratio)。

$$
v=-\frac{\varepsilon_{\mathrm{lateral}}}{\varepsilon_{\mathrm{axial}}}
$$



$$
\varepsilon_{xx}=S_{11}\sigma_{xx}\\\varepsilon_{yy}=S_{21}\sigma_{xx}\\~\\\nu=\frac{\varepsilon_{yy}}{\varepsilon_{xx}}=\frac{S_{21}}{S_{11}}=\frac{-C_{12}}{C_{11}+C_{12}}
$$

## 3 热力学性质

### 3.1 热力学基础

#### 3.1.1 [热力学基本定律](https://en.wikipedia.org/wiki/Laws_of_thermodynamics)

- [**热力学第零定律**](https://en.wikipedia.org/wiki/Zeroth_law_of_thermodynamics)（热平衡定律）
  
  > 如果两个热力学系统中的每一个都与第三个热力学系统处于热平衡 (温度相同)，则它们彼此也必定处于热平衡。
  > 
  > 处于热力学平衡状态的所有物质均具有**某一共同的宏观物理性质**。
  
  热力学第零定律确定了状态函数——*温度*。

- [**热力学第一定律**](https://en.wikipedia.org/wiki/First_law_of_thermodynamics)（能量守恒定律）
  
  > $$
  \Delta U=Q+W
  $$
  > 
  > 物体内能的增加等于物体吸收的热量和对物体所作的功的总和。
  > 
  > [第一类永动机](https://en.wikipedia.org/wiki/Perpetual_motion#Classification)不可能实现。
  
  热力学第一定律确定了状态函数——*内能*和*焓*。

- [**热力学第二定律**](https://en.wikipedia.org/wiki/Second_law_of_thermodynamics)（熵增原理）
  
  > $$
  \eta=\frac A{Q_1}=1-\frac{Q_2}{Q_1}<1
  $$
  > 
  > 克劳修斯表述：热量不能自发地从低温物体转移到高温物体。
  > 
  > 开尔文表述：不可能从单一热源取热使之完全转换为有用的功而不产生其他影响。
  > 
  > 熵增原理：孤立系统的熵永不自动减少，熵在可逆过程中不变，在不可逆过程中增加。
  > 
  > [第二类永动机](https://en.wikipedia.org/wiki/Perpetual_motion#Classification)不可能实现。

- [**热力学第三定律**](https://en.wikipedia.org/wiki/Third_law_of_thermodynamics)（能斯特定理）
  
  > $$
  \lim_{T\to0K}\left(\Delta S\right)_{T}=0
  $$
  > 
  > 普朗克表述：在温度趋于0 $K$时，一切*完美晶体*的熵值趋于一个普遍常量，定义为0。
  > 
  > 能斯特表述：如果等温可逆过程的温度接近 0 $K$，与任何经历该过程的凝聚系统相关的熵变化接近于零。
  > 
  > **不可达原则**：任何过程，无论多么理想化，都不可能在有限的操作次数内将系统的熵降低到其绝对零值。



#### 3.1.2 微观热运动

根据**绝对零度不可达原则**，在任何有限温度下，原子都会具有速度和动能。

- 固体：晶格振动；

- 液体：布朗运动；

- 气体：分子碰撞。

对于晶体来说，晶格振动是典型的热运动，对晶体热学性能起主要贡献，如： 固体热容、热膨胀、熔化、烧结、热传导等。



#### 3.1.3 热胀冷缩

热胀冷缩是大多数物体具有的一种性质，在一般状态下，物体受热以后会膨胀，在受冷的状态下会缩小。

> - 一般情况下，温度升高，分子的动能增加，平均自由程增加，表现为热胀；温度降低时，分子的动能减小，平均自由程减少，表现为冷缩。
> 
> - 水是一个例外，水分子间存在的**氢键具有方向性**，在一定温度范围内，温度下降，水中的氢键数量增加，导致体积随温度下降体积反而增大。



#### 3.1.4 热力学性质

比热容

比热容是**衡量物质容纳热能力**的物理量，它是单位质量或体积的物体升高一定温度所需的能量。

微观解释：通常用微观粒子运动的激烈程度来表征物体温度，影响这些粒子速率改变的因素就能影响物体的比热。

- 固体的粒子间相互作用力比较强，一个原子的运动很容易影响到其他原子的运动，从而导致整体的升温，比热相对较小；

- 液体粒子间相互作用力比固体小，粒子的运动要影响到其它粒子需要更多的能量，因此比热相对较大。



物态

随着温度的升高，材料一般会经历*固态—液态—气态*的物态转变，发生相变过程。



### 3.2 热膨胀系数

[**热膨胀**](https://en.wikipedia.org/wiki/Thermal_expansion)是物质响应温度变化而改变其形状、面积、体积和密度的趋势。材料热胀冷缩的程度可以通过**热膨胀系数**来衡量。

热膨胀系数通常有以下三种定义方法：

- 线膨胀系数：
  
  $$
  \alpha_L=\frac{1}{L}\frac{\mathrm{d}L}{\mathrm{d}T}
  $$

- 面膨胀系数：
  
  $$
  \alpha_A=\frac{1}{A} \frac{\mathrm{d}A}{\mathrm{d}T}
  $$

- 体膨胀系数：
  
  $$
  \alpha_V=\frac{1}{V}\frac{\mathrm{d}V}{\mathrm{d}T}
  $$

对于各向同性材料，

$$
\alpha_A=2\alpha_L\\~\\\alpha_V=3\alpha_L
$$



#### 3.2.1 计算思路

使用*NPT*系综，计算不同温度下$Cu$的体积，曲线的斜率就是$\Delta V/\Delta T$。

#### 3.2.2 输入文件


```python
%%writefile in.thermal_expansion_coefficent

variable T index 300 400 500 600 700 800 900 1000
label T_loop

units metal
boundary p p p
atom_style atomic

# Create lattice
lattice fcc 3.62
region box block 0 8 0 8 0 8
create_box 1 box
create_atoms 1 box

# Set interatomic potential
pair_style eam
pair_coeff 1 1 Cu_u3.eam

# Reset timestep
reset_timestep 0

# Initialize velocity
velocity all create ${T} 87287 dist gaussian

# Equilibrate using NPT ensemble
fix 1 all npt temp ${T} ${T} $(100*dt) iso 0 0 1

# Define thermo output
thermo_style custom step temp epair press lx ly lz
thermo 1000

# Define temperature, length and volume computation
compute actual_T all temp
variable Lx equal lx
variable V equal vol

# Average properties over time and write to a single file
fix 2 all ave/time 100 10 10000 c_actual_T v_Lx v_V file thermal_expansion_data.${T}.txt

# Run simulation
run 10000

# Unfix the NPT ensemble
unfix 1

clear
next T
jump SELF T_loop
```

*拟合数据*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/518c77d1e3d147e597b76d200bd4e3ea/27437a7b-fc82-41ec-8f20-3cdfa1b0cc63.png)

$$
\frac{\Delta L}{\Delta T} = 5.55 \times 10^{-4} \AA / K
$$

### 3.3 比热容

[**比热容**](https://en.wikipedia.org/wiki/Specific_heat_capacity)：单位质量或体积的物质温度升高或降低1 $K$时吸收或放出的热量。

- 体积热容：
  
  $$
  C_{V}=\frac{\Delta E}{\Delta T\cdot V}
  $$

- 质量热容：
  
  $$
  C_{V}=\frac{\Delta E}{\Delta T\cdot V\cdot \rho}
  $$

#### 3.3.1 计算思路

使用*NVT*系综，保持体积在模拟的过程中不变，考虑体系能量随温度的变化。

$$
C_V=\left(\frac{dE}{dT}\right)_V
$$



#### 3.3.2 输入文件


```python
%%writefile in.Cv

units metal 
boundary p p p
atom_style atomic

variable Tstart equal 10.0
variable Tstop equal 1000
variable Tdamp equal 0.2

lattice fcc 3.62
region simbox block 0 8 0 8 0 8
create_box 1 simbox
create_atoms 1 box

pair_style eam
pair_coeff 1 1 Cu_u3.eam

velocity all create ${Tstart} 825577 dist gaussian

# Equilibrium step
fix 1 all nvt temp ${Tstart} ${Tstart} ${Tdamp}
thermo 100
run 5000
unfix 1
reset_timestep 0
thermo 1000
thermo_style custom step temp pe etotal

# Heating step
print "Production steps started!"
fix 2 all nvt temp ${Tstart} ${Tstop} ${Tdamp}
run 120000 
print "Simulation complete!"
```

*拟合数据*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/518c77d1e3d147e597b76d200bd4e3ea/9ce257fd-7e61-4f3e-8138-c4c242a9685d.png)

$$
\frac{\Delta E}{\Delta T} = 0.528 eV/K
$$

### 3.4 相变

晶体从固态转变为液态的过程称为**熔化**，对应的温度为**熔点**；从液态转变为固态的过程称为**凝固**，对应的温度为**凝固点**。

晶体的熔化和凝固是典型的[**一级相变**](https://en.wikipedia.org/wiki/Phase_transition#Modern_classifications)，相变过程中，两相的*化学势连续*，但*化学势的一阶导数不连续*。

$$
\mu_1(T,p)=\mu_2(T,p)\\~\\\frac{\partial\mu_1}{\partial T}\neq\frac{\partial\mu_2}{\partial T},\frac{\partial\mu_1}{\partial p}\neq\frac{\partial\mu_2}{\partial p}
$$

按照$Gibbs$自由能定义，体系的*熵*和*体积*将存在突变，伴随着相变潜热的发生。

#### 3.4.1 计算思路

使用*NPT*系综，通过监控*体积变化*、*焓变*或者其他表示*有序性参量（熵）的值的变化*来标识相变过程，从而得到晶体的熔点或凝固点。

具体而言，主要有以下**转变指标**：

- 体积—温度变化（$V - T$）

- 均方位移（MSD）

- 径向分布函数（RDF）

- 林德曼指数（*Lindemann index*）



#### 3.4.2 熔化与凝固

##### 3.4.2.1 熔化过程


```python
%%writefile in.fusion

units lj
boundary p p p
atom_style atomic

lattice fcc 1.073
region box block 0 10 0 10 0 10
create_box 1 box
create_atoms 1 box
mass 1 1.0

pair_style lj/cut 5.0
pair_coeff 1 1 1.0 1.0 5.0
velocity all create 0.01 87287

timestep 0.005

# Equilibrium step
fix 1 all nvt temp 0.01 0.01 0.2
run 50000
unfix 1

# Heating step
reset_timestep 0
thermo 1000
thermo_style custom step temp pe etotal press vol
fix 1 all npt temp 0.01 0.85 1.0 iso 0 0 1.0
dump 1 all atom 10000 fusion.lammpstrj
run 1000000
```

*熔化过程演示*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/da67be899d33469ba0405ab9411d9c77/ff392bf3-063b-4323-9061-761836485c20.gif)

*V-T图*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/dbe40ddd44a548ecb276094754114c6b/1dd8f53e-c0bd-4412-a97a-ed59ac0eab1d.png)

计算结果表明，分子动力学模拟熔化过程会产生**过热**现象。

##### 3.4.2.2 凝固过程



```python
%%writefile in.quench

units lj
boundary p p p
atom_style atomic

lattice fcc 1.073
region box block 0 10 0 10 0 10
create_box 1 box
create_atoms 1 box
mass 1 1.0

pair_style lj/cut 5.0
pair_coeff 1 1 1.0 1.0 5.0
velocity all create 0.01 87287

timestep 0.005

# Equilibrium step
fix 1 all nvt temp 0.85 0.85 0.2
run 50000
unfix 1

# Heating step
reset_timestep 0
thermo 1000
thermo_style custom step temp pe etotal press vol
fix 1 all npt temp 0.85 0.01 1.0 iso 0 0 1.0
dump 1 all atom 10000 quench.lammpstrj
run 1000000
```

*凝固过程演示*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/da67be899d33469ba0405ab9411d9c77/c1384d7c-7584-4c29-82af-b1b4d8b242f1.gif)

*V-T图*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/da67be899d33469ba0405ab9411d9c77/62a7052e-98d9-4718-9f90-2479164d26ea.png)


计算结果表明，分子动力学模拟凝固过程会产生**过冷**现象。



>  **改进**：在*NPT*系综下进行升温熔化或降温凝固的过程，在特定温度下保持一段时间，以达到热力学平衡，可以部分抑制过热或过冷现象。
> 
> 
> 
> **讨论**：分子动力学方法模拟熔化或凝固过程中，往往会发生*过热*或*过冷*。
> 
> **形核**是理解相变（包括熔化和凝固）的关键概念。以熔化过程为例，形核是指在固体中形成小的液相团簇。这些团簇作为“种子”进行生长，最终导致固体完全熔化。*形核必须克服一个能量势垒。*
> 
> - **临界尺寸**：存在一个临界尺寸，超过该尺寸的核是能量上有利于生长的，而小于该尺寸的核会收缩并消失。*在典型的分子动力学模拟空间尺度形成大于临界尺寸的核是一个小概率事件。*
> 
> - **热波动**：核变通常是由于热波动造成的，这些波动偶尔会在一个小区域内集中足够的能量来克服核变势垒。*在典型的分子动力学模拟时间尺度内发生热波动是一个小概率事件。*
> 
> - 形核位点：在真实材料中，杂质、缺陷或表面通常作为降低核变势垒的形核位点。*这些通常不出现在模拟的理想化条件中。*
> 
> **结论：在分子动力学模拟涉及的时间和空间尺度下，相变初期的形核更加困难。**



#### 3.4.3 体相材料的熔点

##### 3.4.3.1 输入文件


```python
%%writefile in.Cu_fusion

units metal
boundary p p p
atom_style atomic

lattice fcc 3.62
region box block 0 8 0 8 0 8
create_box 1 box
create_atoms 1 box

pair_style eam
pair_coeff 1 1 Cu_u3.eam

thermo 1000
thermo_style custom step temp pe etotal press vol

velocity all create 100.0 4928459 dist gaussian
fix 1 all npt temp 100.0 100.0 0.1 iso 0.0 0.0 1.0
run 10000
unfix 1

reset_timestep 0
compute 1 all msd
fix msd all ave/time 10 100 1000 c_1[1] c_1[2] c_1[3] c_1[4] file msd.dat     # Output msd
fix 1 all npt temp 100 2000 0.1 iso 0 0 1.0
dump 1 all atom 1000 dump.lammpstrj
run 120000
```

##### 3.4.3.2 熔点计算

**均方位移MSD (Mean Square Displacement) ，标识原子偏离其平衡位置的程度**：

$$
MSD=\langle|r(t)-r(0)|^2\rangle 
$$

均方位移MSD与原子的扩散系数存在对应的关系：

- 固态材料的均方位移存在上限；

- 液态均方位移与温度呈线性关系，且其斜率与原子扩散系数存在如下关系：
  
  $$
  {D}=\lim_{t\to\infty}\frac1{2\sigma t}\langle|\boldsymbol{r}(t)-\boldsymbol{r}(0)|^2\rangle 
  $$
  
  其中，$\sigma$为扩散维数。

#### 3.4.4 纳米材料

##### 3.4.4.1 纳米颗粒

纳米颗粒相比于体相结构表面能高，比表面原子数多，表面原子近邻配位不全，活性大于块体材料，因此熔化时所需增加的内能小得多，熔点更低。

- *Cu-Ni纳米颗粒*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/ac1147c726b649a293c2360a1772ee85/03dffa66-9bc1-4f56-9d65-82bf731f1850.png)

- *Pt纳米颗粒*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/811137b4d268442ca4adda97be726b09/72731eaf-9872-47ea-84dd-96e78baee84b.png)

熔点的尺寸效应：

熔点随尺寸减小而降低的现象通常被称为*Melting point depression*。熔点降
低的主要原因是由于纳米颗粒越小，表面原子所占的比重就越大。

$$
T_m(r)=T_m(\infty)\left(1-\alpha\frac dD\right)\\=T_m(\infty)-C/r
$$



相关物理模型：

- Liquid Drop Model

- Liquid Shell Nucleation Model

- Surface Phonon Instability Model

- Bond Order Length Strength Model



**烧结**：颗粒中的原子漫过颗粒边界融合形成一个整体的过程。一些高熔点的材料通常使用烧结作为成型工艺。金属纳米颗粒通常被用作催化剂，其优异催化性能来源于大的比表面积。烧结降低了催化剂的表面积，改变了催化剂的表面结构，是催化活性丧失的重要原因，因此抗烧结是催化剂中是亟需解决的重大问题。





*Ref*: 

[*Nat Commun* **10**, 2583 (2019)](https://doi.org/10.1038/s41467-019-10713-z)

[*Sci Rep* **11**, 19297 (2021)](https://doi.org/10.1038/s41598-021-98704-3)

[*Materials* **2020**, *13*(7), 1507](https://doi.org/10.3390/ma13071507)



##### 3.4.4.2 林德曼指数

[林德曼指数](https://en.wikipedia.org/wiki/Lindemann_index)是原子或分子中热驱动无序的简单量度，一般用于计算**非周期性
体系**的熔点，如纳米颗粒。它定义为键长波动的相对均方根$\delta$：

$$
q_i=\frac{1}{N-1}\sum_{j\neq i}\frac{\sqrt{\langle r_{ij}^2\rangle-\langle r_{ij}\rangle^2}}{\langle r_{ij}\rangle}
$$

或

$$
\delta=\frac2{N(N-1)}\sum_{i<j}\frac{\sqrt{\left\langle r_{ij}^2\right\rangle_t-\left\langle r_{ij}\right\rangle_t^2}}{\left\langle r_{ij}\right\rangle_t}
$$



*Ref*:

[*J. Phys. Chem. A* 2006, 110, 4, 1518–1523](https://doi.org/10.1021/jp053318s "DOI URL")


```python
from ase import io
import numpy as np


def lindemann(filename, n_start=None, n_end=None):
    """
    Calculate the Lindemann index for atomic vibrations using the provided trajectory.

    Parameters:
    - filename: Name of the file containing atomic positions (typically a LAMMPS dump file).
    - n_start: Starting index for the trajectory.
    - n_end: Ending index for the trajectory.

    Returns:
    - Lindemann index value.
    """

    # Construct the slice string based on provided indices
    if n_start is None and n_end is None:
        slice_str = ":"
    elif n_start is None:
        slice_str = f":{n_end}"
    elif n_end is None:
        slice_str = f"{n_start}:"
    else:
        slice_str = f"{n_start}:{n_end}"

    # Read atomic positions from the file for the given range
    atoms = io.read(filename, slice_str, 'lammps-dump-text')

    # Total number of atoms and total number of timesteps
    N = atoms[0].positions.shape[0]
    N_t = len(atoms)
    print('number of atoms in each frame =', N, 'number of frames =', N_t)

    # Initialize arrays to store distances between atoms
    r_square = np.zeros((N_t, N, N))

    # Loop through each timestep
    for t in range(N_t):
        G = np.dot(atoms[t].positions, atoms[t].positions.T)
        H = np.tile(np.diag(G), (N, 1))

        # Compute the squared distance between each pair of atoms at time t
        r_square[t, :, :] = H + H.T - 2 * G

    # Compute average distance and squared average distance over all timesteps
    r_ave = np.mean(np.sqrt(r_square), axis=0)
    print(r_square.shape)
    print(r_ave.shape)
    print('r_ave_max =', np.max(r_ave), 'r_ave_min =', np.min(r_ave))

    r2_ave = np.mean(r_square, axis=0)
    print('r2_ave =', np.max(r2_ave))

    # Extract upper triangular part of matrices
    r_ave_upper = r_ave[np.triu_indices(N, k=1)]
    print('r_ave_upper_max =', np.max(r_ave_upper), 'r_ave_upper_min =', np.min(r_ave_upper))
    r2_ave_upper = r2_ave[np.triu_indices(N, k=1)]

    # Calculate Lindemann criterion for upper triangular elements
    value_upper = np.sqrt(r2_ave_upper - np.power(r_ave_upper, 2)) / r_ave_upper

    # Sum over all unique atom pairs
    total_value = np.sum(value_upper)

    # Normalize by the total number of atom pairs
    return total_value * 2 / (N * (N - 1))
```

##### 3.4.4.3 输入文件

- 纳米颗粒熔化过程


```python
%%writefile in.nanoparticles_fusion

# Initialization
units metal
boundary p p p
atom_style atomic

# Generate NiCu nanoparticles
variable A0 equal 3.589
lattice fcc ${A0}
region mybox block 0 40 0 40 0 40
create_box 2 mybox                         # 2 types of atoms
region CuNi_nano sphere 20 20 20 2
create_atoms 1 region CuNi_nano
set type 1 type/fraction 2 0.5 7777
mass 1 63.54600000                         # Cu
mass 2 58.69340000                         # Ni

# Using EAM potential
pair_style eam/alloy
pair_coeff * * CuNi.eam.alloy Cu Ni
write_data nanoparticle.cif
run 0

# Thermal equilibrium step
thermo 1000
variable j loop 0 20
label loop_j
variable temperature equal 900+10*$j
variable T equal temp
variable Eatom equal etotal/atoms
fix 1 all nvt temp ${temperature} ${temperature} 0.1
run 50000
unfix 1

# Statistical step
fix 2 all ave/time 100 5 1000 v_T v_Eatom file data_ave${temperature}.txt
dump 1 all atom 5000 fusion_${temperature}.atom
fix 1 all nvt temp ${temperature} ${temperature} 0.1
run 500000

unfix 1
undump 1
next j
jump SELF loop_j
```

*纳米颗粒熔化过程演示*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/ac1147c726b649a293c2360a1772ee85/2b01860b-1d1b-4289-ae20-a8b55ec47013.gif)

*Lindemann index - T图*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/ac1147c726b649a293c2360a1772ee85/3179fdd3-dbc4-461d-88ad-686aa040d127.png)


- 纳米颗粒烧结过程


```python
%%writefile in.nanoparticles_sintering

# Initialization
units metal
boundary p p p
atom_style atomic
timestep 0.001

# Get three Pt nanoparticles 
variable A0 equal 3.9239
lattice fcc ${A0}
region mybox block 0 20 0 20 0 20
region sphere_pt1 sphere 7 10 8 3
region sphere_pt2 sphere 10 10 13.4 3
region sphere_pt3 sphere 13 10 8 3
create_box 1 mybox
create_atoms 1 region sphere_pt1
create_atoms 1 region sphere_pt2
create_atoms 1 region sphere_pt3

# Using EAM potential
pair_style eam
pair_coeff 1 1 Pt_u3.eam

# Output initial structure
dump 1 all cfg 1 coord.*.cfg mass type xs ys zs
run 0
undump 1

# Define triple_neck region
region triple_neck block 7.5 12.5 5 15 9 13
group 1 dynamic all region triple_neck every 1000

# Output number of atoms in triple neck region
variable N equal step
variable T equal temp
variable Natom equal count(all)
variable V equal vol/v_Natom
variable sinter_atom equal count(1)
dump 1 all xyz 1000 melt.xyz
thermo 1000
fix extra all print 1000 "${N} ${T} ${sinter_atom}" file data.txt

# Run in 500k, 1000k, 1400k, respectively
fix 1 all npt temp 500 500 0.1 iso 1 1 1
run 10000
unfix 1
fix 1 all npt temp 1000 1000 0.1 iso 1 1 1
run 20000
unfix 1
fix 1 all npt temp 1400 1400 0.1 iso 1 1 1
run 30000
unfix 1
```

*纳米颗粒烧结过程演示*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/811137b4d268442ca4adda97be726b09/2ae29c08-489b-43ef-abce-28d1bb580d63.gif)

*Temp/Atom_number - Steps图*：

![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/811137b4d268442ca4adda97be726b09/d134ebfa-ede2-406a-98af-cb7ca87b55c2.png)

## 4 输运性质

### 4.1 离子输运（能源材料）

通过模拟原子、分子等粒子在其相互作用力下的运动过程，探索材料结构与动力学性质，可广泛应用于能源材料的原子尺度相互作用机制，包括结合能规律、变形机制、热输运规律、偏析相形成等与能源材料密切相关的行为，在能源材料研究中具有重要意义：

- 材料相的形成及性能预测：如相转变、力学强度、热导率、体积膨胀等。

- 界面和表面研究：如太阳能电池、储能设备和电化学电池。

- 离子传输和扩散：如锂离子电池、燃料电池和电解质材料。

- 储氢材料研究：可模拟氢气在储氢材料中的吸附、扩散和释放过程。

#### 4.1.1 Li-S电池膨胀

**锂离子电极电池基本组成**：

- 阳极（负极）
  
  放电时由外电路*失去*电子；发生**氧化反应**；电位*低*。

- 隔膜

- 阴极
  
  放电时由外电路*获得*电子；发生**还原反应**；电位*高*。

- 有机电解液

- 电池外壳

锂电池经过高温存储和循环，容易发生鼓胀，由此产生一系列安全问题：

1. 由于电极材料膨胀可能导致电池部件的*物理损坏*；

2. 体积膨胀可能导致正极和负极之间的物理接触，*引发电池短路*；

3. 可能导致电极材料的失效，*降低电池的容量和循环寿命*。

**LiS电池**

- 优点
  
  - 高比容量
  
  - 高能量密度
  
  - 低成本
  
  - 环境友好

- 缺点
  
  - 电导率低，放电效率低
  
  - 多硫离子溶于电解液，循环性能差
  
  - 体积膨胀大，循环寿命短

*Ref*:

- [*Acc. Chem. Res.* 2013, 46, 5, 1135–1143](https://pubs.acs.org/doi/10.1021/ar3001348)

- [*J. Phys. Chem. Lett.* 2013, 4, 19, 3227–3232](https://pubs.acs.org/doi/10.1021/jz401763d)

##### 4.1.1.1 模拟思路

1. 随机产生一定数目的$Li$原子坐标填入$S$结构中，作为$$Li_xS_{128}$初始结构；

2. $Li$嵌入后$S$变成无定形相，考虑**升温熔化+快速退火**构建无定形结构；

3. 在*NPT*系综和标准大气压下统计室温的$Li_xS$无定形相的体积，与$Li$嵌入前的体积对比计算体积膨胀率。

**势函数**：[***Phys. Chem. Chem. Phys.***, 2015,**17**, 3383-3393](https://pubs.rsc.org/en/content/articlelanding/2015/cp/c4cp04532g)

##### 4.1.1.2 输入文件


```python
%%writefile in.thermal_expansion

units real
boundary p p p
atom_style charge

read_data Li2S.dat

pair_style reax/c NULL checkqeq no
pair_coeff * * ffield.reax.LiS S Li

# Fix 1 all qeq/reax 1 0.0 10.0 1e-6 reax/c
neighbor 3.0 bin
min_style cg
minimize 1e-4 1e-6 100 1000

# Stablize at room Temp
velocity all create 300 96588
thermo 1000
thermo_style custom step temp press density vol
compute myRDF all rdf 50 1 1 1 2 2 1 2 2
fix 2 all ave/time 100 10 10000 c_myRDF[*] file initial.rdf mode vector
fix 3 all nvt temp 300 300 10
run 10000
unfix 2
unfix 3

# Heating slowly
fix 4 all npt temp 300.0 1800.0 100 iso 1.0 1.0 1000 drag 1.0
run 50000
unfix 4

# Keep melting
fix 5 all npt temp 1800.0 1800.0 100 iso 1.0 1.0 1000 drag 1.0
fix 6 all ave/time 100 10 5000 c_myRDF[*] file melted.rdf mode vector
run 5000
unfix 5
unfix 6

# Rapidly quenching
fix 7 all npt temp 1800.0 300.0 1 iso 1.0 1.0 1000 drag 0.1
thermo 10
run 100
unfix 7

# Stablize the amorphous LixS at RT
fix 8 all ave/time 100 10 5000 c_myRDF[*] file amorphous.rdf mode vector
fix 9 all npt temp 300 300 100 iso 1.0 1.0 1000 drag 1.0
thermo 1000

```

> [径向分布函数(Radial Distribution Function, RDF)](https://en.wikipedia.org/wiki/Radial_distribution_function)
> 
> ![alt](https://bohrium.oss-cn-zhangjiakou.aliyuncs.com/article/16392/c7aa6c89e33141a384510c426c263a08/cf983c01-462b-4a03-afcd-11b5d2a054df.png)
> 
> 
> 原子相对密度随半径变化的函数
> 
> $$
g(r)=\frac{\left\langle N(r\pm\frac{\Delta r}2)\right\rangle}{\Omega(r\pm\frac{\Delta r}2)}\cdot\frac1\rho 
$$
> 
> 该物理量可以用来表征原子的**局部结构有序度**，标识有序-无序相的变化。
> 
> - 固体：近邻原子出现在特征距离处；
> 
> - 液体：近邻原子在特征距离附近平缓连续变化；
> 
> - 气体：近邻原子不具备长程和短程有序性。

#### 4.1.2 Mg储氢机制

**储氢现状**：

- 高压缸瓶：密度高，压力大，存在安全隐患；

- 多孔吸附：比重低，吸附弱；

- 化合物氢：比重高，吸附强，难以脱附释放。

##### 4.1.2.1 模拟思路

研究对象：金属Mg(100)面

研究目的：

1. 探究Mg的储氢机制；

2. 探究不同温度对Mg储氢效果的影响。

**势函数**：[***Computational Materials Science*** 154 (2018) 295–302](https://www.sciencedirect.com/science/article/pii/S0927025618304865)

##### 4.1.2.2 输入文件


```python
%%writefile in.hydrogen_storage

units metal
dimension 3
atom_style atomic
boundary p p p
read_data MgH.data
neighbor 0.5 bin
neigh_modify every 1 delay 0 check yes

region Mg block INF INF INF INF -0.01 5 units box
group Mg region Mg
velocity Mg set 0 0 0
fix 01 Mg setforce 0 0 0
region H2 block INF INF INF INF 6 INF units box
group H2 region H2

mass 1 24.3050
mass 2 1

pair_style adp
pair_coeff * * Mg_H.adp.alloy.txt Mg H

region 1H block INF INF INF 9 INF units box
group mobile1 dynamic all region 1H every 100
variable number1 equal count(mobile1)

region 2H block INF INF INF 5 9 units box
group mobile2 dynamic all region 2H every 100
variable number2 equal count(mobile2)

region 3H block INF INF INF -0.001 5 units box
group mobile3 dynamic all region 3H every 100
variable number3 equal count(mobile3)

# compute 1 Mg temp
compute 2 H2 temp
timestep 0.001

variable T equal TEMP
thermo 1000
thermo_style custom step temp pe etotal press vol c_2

min_style cg
minimize 1e-5 1e-5 10000 10000
run 0

velocity H2 create ${T} 71717
dump 1 all atom 100 dump.xyz
fix 1 H2 nvt temp ${T} ${T} 0.1
thermo_modify lost ignore
fix extra all print 100 "${N} ${P} ${number1} ${number2} ${number3}" file energy.data
run 200000
```

### 4.2 热输运（热电材料）

#### 4.2.1 热传导理论基础

热传递机制：

- 固体：热传导

- 液体：热对流

- 气体：对流和辐射

热传导：热量从固体材料高温部分向低温部分转移的过程。

根据热力学第二定律（克劳修斯表述），在自然条件下热量只能从高温物体向低温物体转移，这个转变过程是**自发**，**不可逆**的。

由热力学第三定律，原子始终在格点平衡位置附近进行热振动。**材料各种热学性能的物理本质，均与其晶格热振动有关。**

绝缘材料中的导热现象主要由原子（晶格）的振动引起

在量子力学中：

$$
H=\frac12m\dot{x}^2+\frac12kx^2=\frac12\dot{q}^2+\frac12\omega^2q^2
$$

谐振子具有分立的本征能量：

$$
\varepsilon_n=\left(n+\frac12\right)\hbar\omega,\quad n=0,1,\cdots,\infty 
$$

定义[声子](https://en.wikipedia.org/wiki/Phonon)为晶格振动的能量量子，则声子能量$\propto\hbar\omega=\hbar\sqrt{\frac km}$，声子之间的碰撞导致热量的传递，因此原子/分子间作用力（键）越强，导热能力越强。

声子碰撞过程：

- Normal process（正常过程）：无热阻

$$
\vec{q}_1+\vec{q}_2=\vec{q}_3
$$

- Umklap process（倒逆过程）：有热阻

$$
\vec{q}_1+\vec{q}_2=\vec{q}_3+\vec{K}_h
$$

> - 长波声子具有更长的波长，易发生倒逆过程，因此对热传导性质影响更大。
> 
> - 声学声子对热导率的贡献远大于光学声子。



#### 4.2.2 热导率

[热导率](https://en.wikipedia.org/wiki/Thermal_conductivity_and_resistivity)，又称导热系数，是物质导热能力的量度。

标准单位：$W/(m\cdot K)$

$$
\kappa=\frac13C\bar{v}\bar{\lambda}
$$

- 与**固体比热**$C$成正比。因为比热越大，每个固体分子可携带的热量也越多；

- 与**声子平均速度**$\bar{v}$成正比。因为声子速度越快，热传递也越快；

- 与**声子平均自由程**$\bar{\lambda}$（与声子碰撞频率相关）成正比。$\bar{\lambda}$越大，声子携带热量传递越远。

上式中，声子平均速度与温度变化的关联性不大，多数情况可以作为常数处理；平均自由程较难计算，与温度和声子碰撞的方式、频率等因素有关。

**晶体热容**：

$$
C_V=\left(\frac{\partial U}{\partial T}\right)_V
$$

固体热容主要来源于两部分：

- 晶格振动（声子运动）——晶格热容$C_V^l$

- 电子热运动——电子热容$C_V^e$

$$
C_V=C_V^l+C_V^e
$$

通常情况下，$C_V^e<<C_V^l$，在极低温时，则需要考虑电子热容贡献。

***实验规律***：

1. 高温时，晶体热容
   
   $$
   C_V=3Nk_B=3\nu R=24.9J/K\cdot mol
   $$

2. 低温时，
   
   绝缘体热容$\propto T^3$；
   
   导体热容$\propto T$。

晶格振动的能量是量子化的，频率为$\omega$的声子能量为

$$
E_n=(n+\frac12)\hbar\omega
$$

其中，$\frac12\hbar\omega$代表零点振动能，对热容没有贡献，因此

$$
E_n=n\hbar\omega
$$

其中$n$是频率为$\omega$的谐振子的平均声子数，根据**玻色统计理论**，

$$
n(\omega)=\frac1{e^{\hbar\omega/k_BT}-1}
$$

因此，温度为$T$时，频率为$\omega$的振动能量

$$
\bar{E}(\omega)=\frac{\hbar\omega}{e^{\hbar\omega/k_BT}-1}
$$

晶体由$N$个原子组成，每个原子有$3$个自由度，共有$3N$个分立的振动频率，晶体内能

$$
U=\sum_{i=1}^{3N}\bar{E}\left(\omega_i\right)=\sum_{i=1}^{3N}\frac{\hbar\omega_i}{e^{\hbar\omega_i/k_BT}-1}
$$

- **爱因斯坦模型**假设所有原子以同样的频率振动，

$$
U=\sum_{i=1}^{3N}\bar{E}\left(\omega_i\right)=\frac{3N\hbar\omega}{e^{\hbar\omega/k_BT}-1}
$$

该模型在高温下比热为$3Nk_B$，与经典结论一致，但低温下偏差较大。

- **德拜模型**将频率分布用积分函数表示，

$$
\int_0^{\omega_m}g(\omega)d\omega=3N
$$

因此，

$$
U=\sum_{i=1}^{3N}\frac{\hbar\omega_i}{e^{\hbar\omega_i/k_BT}-1} =\int_0^{\omega_m}\frac{\hbar\omega}{e^{\hbar\omega/k_bT}-1}g(\omega)d\omega 
$$

热容

$$
C_V=\left(\frac{\partial U}{\partial T}\right)_V=\int_0^{\omega_m}k_B\left(\frac{\hbar\omega}{k_BT}\right)^2\frac{e^{\hbar\omega/k_BT}g(\omega)d\omega}{(e^{\hbar\omega/k_BT}-1)^2}
$$

变量代换

$$
f_D\left(\frac{\theta_D}T\right)=3\left(\frac T{\theta_D}\right)^3\int_0^{\frac{\theta_D}T}\frac{e^xx^4}{(e^x-1)^2}dx
$$

则低温下，$C_V{\sim}\frac{12\pi^4R}5{(\frac T{\theta_D})}^3$

**声子平均自由程**与声子数密度成反比，

$$
\bar{\lambda}\propto\frac1n\propto e^{\frac\alpha T}
$$

- 低温下，声子平均自由程$L$近似为粒径大小，$\kappa\propto T^3$；

- 随着温度增加，热容增加，声子的平均自由程减小，二者竞争使得$\kappa$出现极大值；

- 高温下，声子平均自由程$L$随温度升高减小，声子热容$C_V$趋于常数，$\kappa$下降。

**热电材料**

- [**热电效应**](https://en.wikipedia.org/wiki/Thermoelectric_effect)
  
  - [**塞贝克（Seebeck）效应**](https://en.wikipedia.org/wiki/Thermoelectric_effect#Seebeck_effect)
  
  - [**珀尔贴（Peltier）效应**](https://en.wikipedia.org/wiki/Thermoelectric_effect#Peltier_effect)

- 热电优值

$$
ZT=\frac{S^2T\sigma}K
$$

[***Fourier's law***](https://en.wikipedia.org/wiki/Thermal_conduction#Fourier's_law)

> 类比[Ohm's law](https://en.wikipedia.org/wiki/Ohm%27s_law)
> 
> **欧姆定律**描述了电子在导电材料内部的输运规律：
> 
> $$
\frac{dq}{dt}=\sigma\cdot A\frac Vl
$$
> 
> **傅里叶定律**描述了热量在材料内部的输运规律：
> 
> $$
\frac{dQ}{dt}=\kappa\cdot A\frac{\Delta T}l
$$



##### 4.2.2.1 计算思路

傅里叶定律：

$$
\frac{dQ}{dt}=\kappa\cdot A\frac{\Delta T}l
$$

计算出单位截面积的热流和温度梯度，求得材料的热导率。



- 采用特定方法对材料某一长度方向施加热流，计算热功率；

- 经过一定时间后材料内部温度达到稳态，统计温度梯度；

- 结合截面积和长度参数，计算热导率。



##### 4.2.2.2 输入文件

- ***固定热流交换法***



```python
%%writefile in.langevin

units lj
atom_style atomic
lattice fcc 0.6
region box block 0 10 0 10 0 20
create_box 1 box
create_atoms 1 box
mass 1 1.0
velocity all create 1.35 71717
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0
neighbor 0.3 bin

region hot block INF INF INF INF 0 1
region cold block INF INF INF INF 10 11

compute Thot all temp/region hot
compute Tcold all temp/region cold

# 1st Equilibration run
fix 1 all nvt temp 1.35 1.35 0.5
thermo 100
run 1000

velocity all scale 1.35
unfix 1

# 2nd Equilibration run
fix 1 all nve
fix hot all heat 1 100.0 region hot
fix cold all heat 1 -100.0 region cold
thermo_style custom step temp c_Thot c_Tcold
thermo 1000
run 10000

compute ke all ke/atom
variable temp atom c_ke/1.5

compute layers all chunk/atom bin/1d z lower 0.05 units reduced
fix 2 all ave/chunk 10 100 1000 layers v_temp file profile.heat
variable tdiff equal f_2[11][3]-f_2[1][3]

fix ave all ave/time 1 1 1000 v_tdiff ave running start 13000
thermo_style custom step temp c_Thot c_Tcold v_tdiff f_ave
run 20000
```

- ***Muller-Plathe（速度交换）法***


```python
%%writefile in.mp

units lj
atom_style atomic
lattice fcc 0.6
region box block 0 10 0 10 0 20
create_box 1 box
create_atoms 1 box
mass 1 1.0
velocity all create 1.35 71717
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0
neighbor 0.3 bin

# 1st Equilibration run
fix 1 all nvt temp 1.35 1.35 0.5
thermo 100
run 1000
velocity all scale 1.35
unfix 1

# 2nd Equilibration run
compute ke all ke/atom
variable temp atom c_ke/1.5
fix 1 all nve
compute layers all chunk/atom bin/1d z lower 0.05 units reduced
fix 2 all ave/chunk 10 100 1000 layers v_temp file profile.mp
fix 3 all thermal/conductivity 10 z 20
variable tdiff equal f_2[11][3]-f_2[1][3]
thermo_style custom step temp epair etotal f_3 v_tdiff
thermo 1000
run 20000

# Thermal conductivity calculation
fix 3 all thermal/conductivity 10 z 20
fix ave all ave/time 1 1 1000 v_tdiff ave running
thermo_style custom step temp epair etotal f_3 v_tdiff f_ave
run 20000
```

##### 4.2.2.3 注意事项

- 利用傅里叶定律计算热导率的关键和难点是**在系统中产生稳定、线性的温度梯度分布**；

- 无法获得温度梯度时，需要**调整系综、热流交换速率和时间步等参数**；

- 温度梯度不宜太大，一般需要**保持高温区和低温区的温差在100K以内**，平均值在给定的系综温度附近，不要有太大漂移。



#### 4.2.3 总结

##### 4.2.3.1 非平衡分子动力学（NEMD）

通过傅里叶定律计算热导率的方法：

- 模拟过程中对象各部分能量处于**非平衡**状态
  
  - 固定热流交换法与外界有人为的能量交换
  
  - MP（速度交换）法各组原子之间有人为的能量交换

- 经过一段时间后，扩散效率与人为交换能量相当时，体系各组的温度趋向稳定

- 基于傅里叶定律计算热导率的方法通常被称为***非平衡分子动力模拟(NEMD)***。



此外，各向异性材料的热导率是一个张量，***NEMD***一次模拟只能获得一个方向的热导率，适合*低维材料*热导率的模拟。



##### 4.2.3.2 平衡态分子动力学（EMD）

基于[Green-Kubo 线性响应理论](https://en.wikipedia.org/wiki/Green%E2%80%93Kubo_relations)，材料的热导率与平衡态下涨落-耗散的热流相关，可以通过平衡态的热流自相关函数来计算热导率。

$$
\kappa_{\mu\nu}(\tau_m)=\frac1{\Omega_{k_B}T^2}\int_0^{\tau_m}\left\langle\bar{J}_\mu(\tau)J_\nu(0)\right\rangle d\tau 
$$

离散化时间步求平均，

$$
\kappa_{\mu v}(\tau_M)=\frac{\Delta t}{\Omega kBT^2}\sum_{m=1}^M(N-m)^{-1}\sum_{m=1}^{N-m}J_\mu(m+n)J_\nu(n)
$$

此即平衡态下分子动力学方法计算热导率的公式。



*ref*: [Phys. Rev. B **65**, 144306](https://doi.org/10.1103/PhysRevB.65.144306)