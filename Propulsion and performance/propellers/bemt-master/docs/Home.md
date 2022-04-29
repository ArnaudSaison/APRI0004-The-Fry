Welcome to the BEMT wiki! This wiki contains all the documentation regarding
this BEMT code.

## Introduction

**BEMT** is an implementation of the Blade Element Momentum Theory in Matlab
developed mainly for teaching purposes at the [University of
Liege](https://www.am.uliege.be/) (Belgium).

The aim of this code is to allow students to quickly test different
rotors/propellers configurations and see how they perform in terms of thrust,
power, torque and efficiency. This **BEMT** code should therefore be considered
primarily as an _analysis_ code. Although, minor modifications could easily
extent it to be an effective design tool[^1].

Although the code is meant for teaching purposes and not as a proper design tool
for industrial purposes, it still comes with lots of [interesting
  features](docs/FeaturesOverview.md).

### Solvers

As its name details it, the Blade Element Momentum Theory stands on two main
group of equations: the blade element equations and the conservation of
momentum. These equations come in different forms: some are linearised, some
involve the definition of induction factors or of specific velocity triangles
and some try to combine everything in a single main equation.

The present code includes three different solvers, each based on a different
formulation of the base equations:

- [Leishman](docs/Solver-Leishman.md)
- [IndFact](docs/Solver-IndFact.md)
- [Stahlut](docs/Solver-Stahlut.md)

### Contents

_If you are browsing this documentation on Gitlab's wiki, you can use the
navigation menu on the right for the different sections._

- Introduction
  - [Features overview](docs/FeaturesOverview.md)
  - [Prerequisites](docs/Prerequisites.md)
  - [Usage](docs/Usage.md)
  - [Configuration](docs/Code-Config.md)
- Solvers
  - [Leishman](docs/Solver-Leishman.md)
  - [IndFact](docs/Solver-IndFact.md)
  - [Stahlut](docs/Solver-Stahlut.md)
- Code architecture
  - [Configuration](docs/Code-Config.md)
  - [Airfoil data](docs/Code-AirfoilData.md)
  - [Pre-processing](docs/Code-Preproc.md)
  - [Post-processing](docs/Code-Postproc.md)
  - [Analysis](docs/Code-Analysis.md)
- Misc
  - [Validation](Validation.md)
  - [Roadmap](Roadmap.md)
  - [References](References.md)
  - [Troubleshooting](Troubleshooting.md)

[^1]: This was purposely left out in order to let the students do the design
  work on their own. If you are interested by a design version, reach out to me
  by email: tlambert@uliege.be.
