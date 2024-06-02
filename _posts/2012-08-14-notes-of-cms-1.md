---
title: 'Notes of CMS 1'
date: 2024-03-10
tags:
  - notes
  - Computational Materials Science
---

《计算材料学》学习笔记
===============

<div style="color:black; background-color:#FFF3E9; border: 1px solid #FFE0C3; border-radius: 10px; margin-bottom:0rem">
    <p style="margin:1rem; padding-left: 1rem; line-height: 2.5;">
        ©️ <b><i>Copyright 2024 @ Authors</i></b><br/>
        <i>作者：
            <b>
            <a href="mailto:your_address@email.com">王江海 📨 </a>
            </b>
        </i>
        <br/>
        <i>日期：2024-03</i><br/>
        <i>共享协议：</a>本作品采用<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">知识共享署名-非商业性使用-相同方式共享 4.0 国际许可协议</a>进行许可。</i><br/>
        <i>快速开始：点击上方的</i> <span style="background-color:rgb(85, 91, 228); color:white; padding: 3px; border-radius: 5px;box-shadow: 2px 2px 3px rgba(0, 0, 0, 0.3); font-size:0.75rem;">开始连接</span> <i>按钮，选择 <b><u>bohrium-notebook:2023-04-07</u>镜像</b> 和任意配置机型即可开始。
    </p>
</div>



## 计算方法分类

### 纳观尺度

**基本粒子**：电子

**主要理论**：量子力学

$$
i\hbar\frac{\partial\psi}{\partial t}=-\frac{\hbar^2}{2m}\frac{\partial^2\psi}{\partial x^2}
$$

**主要方法**：基于量子化学、密度泛函理论（DFT）的第一性原理计算

**常用软件**：VASP、Gaussian、QuantumEspresso

**计算内容**：

- 结构性质：分子、晶体结构预测

- 表面性质：重构、缺陷、表面能

- 力学性质：弹性常数、杨氏模量

- 磁学性质：磁性、自旋轨道耦合

- 电子性质：分子轨道、能带结构、价态

- 光学性质：吸收光谱、折射率

- 动力学模拟：扩散系数、反应动力学过程

### 微观尺度

> 当原子**德布罗意波长**远小于晶格常数时，可以用牛顿经典力学近似描述原子的运动。

**基本粒子**：原子

**主要理论**：牛顿力学

$$
\overrightarrow{F}=m\overrightarrow{a}
$$

**主要方法**：蒙特卡洛（MC）、分子动力学（MD）

**常用软件**：LAMMPS、GROMACS

**计算内容**：

- 力学性质：微观塑性变形机制、拉伸断裂机理

- 热学性质：相变、热膨胀系数

- 物质结构/相互作用：吸附、扩散、缺陷运动、超分子自组装、表面能

### 介观尺度

**基本粒子**：粗粒子（或组织、结构）

**主要理论**：牛顿力学、统计热力学

**主要方法**：

- 格子玻尔兹曼（LBM）

- 耗散粒子动力学（DPD）

- 布朗动力学（BD）

- 位错动力学（DDD）

- 相变动力学（PTD）

**常用软件**：LAMMPS、GROMACS

### 宏观尺度

**基本粒子**：连续体

**主要理论**：理论力学、流体力学

**主要方法**：

- 有限元法（FEM）

- 有限差分法（FDM）

- 有限体积法（FVM）

常用软件：ANSYS、ABAQUS

**计算内容**：

- 冲击损伤仿真

- 力/热/电/磁场耦合分析

- 温度分布模拟

- 切削表面形成过程

- 裂纹的扩展及应力分布等

---

## 计算材料基本流程

1. 确定物理模型（误差来源）

2. 选择数值方法（误差来源）

3. 分析计算结果（数值结果的选取）

4. 得到物理结论