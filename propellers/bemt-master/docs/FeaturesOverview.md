- [Main current features](docs/FeaturesOverview.md#main-current-features)
  - [Geometry](docs/FeaturesOverview.md#geometry)
  - [Solvers](docs/FeaturesOverview.md#solvers)
  - [Aerodynamics](docs/FeaturesOverview.md#aerodynamics)
  - [Roadmap](docs/Roadmap.md)
- [Analysis](docs/FeaturesOverview.md#analysis)
  - [Rotor performance](docs/FeaturesOverview.md#rotor-performance)
  - [Plots](docs/FeaturesOverview.md#plots)
  - [3D-Views](docs/FeaturesOverview.md#3d-views)
  
---

## Main current features

### Geometry

- Linear twist/taper
- Ideal twist/taper (from Leishman definition)
- Custom twist/taper

### Solvers

- Linearised equations for hovering rotors ([Leishman](docs/Solver-Leishman.md))
- Full system with induction factors ([IndFact](docs/Solver-IndFact.md))
- Full system with one single equation ([Stahlut](docs/Solver-Stahlut.md))

### Aerodynamics

#### Losses

- Prandtl model for tip and hub losses

### Roadmap

Lots of future improvements have been thought of. See the full list on the
[Roadmap](docs/Roadmap.md)

## Analysis

After each run, the code outputs the following information by default (this can
be toggled off in the [configuration](docs/Code-Config.md)).

### Performance

- Rotor thrust, torque, power and propulsive efficiency.

### Plots

- Spanwize evolution of local angle of attack
- Spanwize evolution of $`C_T`$, $`C_P`$, $`C_Q`$
- Spanwize evolution of $`T`$, $`P`$, $`Q`$
- Spanwize evolution of $`C_l`$, $`C_d`$
- Spanwize evolution of $`L`$, $`D`$

### 3D Views

- A view of a single blade, colored to show the contribution to the total thrust
- A view of the complete rotor/propeller, colored to show the contribution to
  the total thrust
