FemLab
======

FemLab is a Finite Element library written in Julia language.
The purpose of this library is to aid the research of new algorithms for the finite element method.

Currently this library solves static and dynamic analysis in two and three-dimensions. It is also capable to simulate cracking based on the
discrete crack approach. There is work on the way to implement a dynamic solver

**Available elements:**
Truss (2D and 3D), Beam (2D), Solid elements (2D and 3D). Also includes embedded and semi-embedded rod
elements used for the simulation of reinforced media.

**Constitutive models:**
 Elastic, Perfectly Plastic, Mazars Damage, Drucker Prager, Plastic joint, CEB joint, etc.

Get FemLab:
Pkg.clone("https://github.com/RaulDurand/FemLab.jl")
