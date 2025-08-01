---
title: 'ATOMS - Installation & Use'
date: 2025-2-28
tags:

  - ATOMS
---

# ATOMS - Installation & Use

## ATOMS 1: Download, Installation & Testing

1. [Shape Software](https://www.shapesoftware.com/)
2. Download Demos $\Rightarrow$ ATOMS for Windows Demo and download
3. Install ATOMS in preferred folder
4. Cancel registration and use in Demonstration Mode
5. Open a test file (\*.str) from the 'installation_folder'/samples/inorgan/\*.str

## ATOMS 2: Manual Structure Input

Input data to draw the rutile structure

1. Create a New Data set $\Rightarrow$ Add Structure axes, Cell parameters, Title

2. Add Space group from table $\Rightarrow$ Highlight + Select + Apply symbol

3. Set plot Boundary Option to Default Unit Cell $\Rightarrow$ Select -1 to 1 inclusive

4. Select Crystal Forms for Display

5. Add Input Atoms (Name, Type, Radius, x, y, z) $\Rightarrow$ Select Fill Color

6. Polyhedra, Bonds and Display Crystal Edges and Faces

   Perspective viewing, Stereopair and Rims

   Linewidths, Shading and Background Color

   Crystal Axes Display, Unit Cell Display and Thermal Ellipsoid Parameters

   Initial Orientation, Scaling and Centering

7. Calculation Output, Save as and Calculate now

## ATOMS 3: Atom Color, Scaling & Atom Shading

1. Change atom color

   Input 1 $\rightarrow$ Atoms $\rightarrow$ Revise $\rightarrow$ Select Fill Color

2. Change background color

   Input 2 $\rightarrow$ Background color $\rightarrow$ Select Color

3. Scaling

   Input 2 $\rightarrow$ Scaling $\rightarrow$ Fixed Scale Factors $\rightarrow$ Screen (typically 0.3 to 0.6)

   <u>or</u> use Rescale on right-hand tool bar

4. Shading

   Input 2 $\rightarrow$ Shading $\rightarrow$ Check/Uncheck Shading for atoms/bonds/polyhedra $\rightarrow$ Number of Zones (usually 32)

## ATOMS 4: Crystallographic Orientation & Unit Cell

1. Insert Crystal Axes

   Input 2 $\rightarrow$ Crystal Axes $\rightarrow$ Crystal Axes Display $\rightarrow$ Show Axes

   (Crystal Axes Display) Font $\rightarrow$ Select Style and Size

2. Show Unit Cell

   Input 2 $\rightarrow$ Unit Cell $\rightarrow$ Show Unit Cell

3. Specify Orientation

   Rotation $\rightarrow$ Align Face or Vector $\rightarrow$ Vector select $\rightarrow$ Use default alignment face

4. Manual Rotation

   Use arrows left-hand toolbox

5. Check Alignment

   Rotation $\rightarrow$ Current Orientation

## ATOMS 5: Crystallographic, Perspective & Clinographic Projection

1. Crystallographic Projection

   Input 2 $\rightarrow$ Perspective $\rightarrow$ Perspective viewing <u>unchecked</u>

2. Perspective Viewing

   Input 2 $\rightarrow$ Perspective $\rightarrow$ Perspective viewing <u>checked</u>

3. Clinographic Projection

   Check box in right-hand panel

## ATOMS 6: Extracting Bond Lengths & Drawing Polyhedra

1. Geometric Projection

   Input 1 $\rightarrow$ Atoms $\rightarrow$ Coordination $\rightarrow$ Coordination of Atoms

   Set Distance limit, Central-ligand (typically 2 - 3 $\AA$)

2. Analyze Bond List

   Note Central atom type number + bonding atom type number + number of bonding atoms + longest bond length

3. Polyhedra Specification

   Input 1 $\rightarrow$ Polyhedra $\rightarrow$ Polyhedra Data

   Add polyhedron number 1 $\rightarrow$ Coordination number $\rightarrow$ Maximum bond distance

   Types Central $\rightarrow$ Type Ligand

   Colors $\rightarrow$ Select Fill Color $\rightarrow$ OK $\rightarrow$ Calculate (left panel)

## ATOMS 7: Boundary Specification

1. Default Boundary Settings

   Input 1 $\rightarrow$ Boundary $\rightarrow$ Unit Cell $\rightarrow$ Check

2. Manual Boundary Settings using Crystallographic Forms

   Input 1 $\rightarrow$ Boundary $\rightarrow$ Unit Cell $\rightarrow$ Change Option $\rightarrow$ Boundary Option $\rightarrow$ Enter Forms $\rightarrow$ OK $\rightarrow$ Revise Forms $\rightarrow$ Central Distance $\rightarrow$ Adjust distance in direction to expand  or reduce plotting

3. Set Drawing Center

   Right-hand toolbox $\rightarrow$ Set Center $\rightarrow$ Define center with mouse

## ATOMS 8: Drawing Bonds & Translucent Polyhedra

1. Determine Bonding Atoms and Bond Lengths

   Input 1 $\rightarrow$ Atoms $\rightarrow$ Input Atoms $\rightarrow$ Coordination of Atoms $\rightarrow$ Distance limit, central-ligand (2-3) $\rightarrow$ OK $\rightarrow$ Inspect distances for bond lengths

2. Specify Bonds and Draw

   Input 1 $\rightarrow$ Bonds $\rightarrow$ Add Bonds $\rightarrow$ Insert Atom type

   Specify maximum bond distance envelop $\rightarrow$ Specify bond fill color $\rightarrow$ OK $\rightarrow$ Replot

3. Define Polyhedra

   Input 1 $\rightarrow$ Polyhedra $\rightarrow$ Add Polyhedra $\rightarrow$ Coordination number + Maximum bond distance

   Specify Central atom and Ligand $\rightarrow$ Specify Fill Color

   Input 2 $\rightarrow$ Check shading for polyhedra

4. Define Translucent Polyhedra and Bonds

   Draw Mode (left-hand tool bar) $\rightarrow$ OGL Single

   Input 1 $\rightarrow$ Polyhedra $\rightarrow$ 3D Parameters

   Specify Faces as Translucent (Opacity 0.6) $\rightarrow$ Specify Central-Ligand Bonds as Rods (0.05) $\rightarrow$ OK

## ATOMS 9: Crystallographic Information File (CIF) Upload & Use

1. File $\rightarrow$ Import File $\rightarrow$ Select CIF
2. Boundary Option $\rightarrow$ Default Unit Cell $\rightarrow$ Open \*.cif $\rightarrow$ Save as \*.str $\rightarrow$ Calculate now $\rightarrow$ Yes
3. Set atom colors and size to preference
4. Draw unit cell, adjust boundary and resize
5. Check bond distances and coordination and draw polyhedra
6. Adjust orientation and insert crystallographic axes
