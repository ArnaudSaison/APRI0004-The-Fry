- [Disclaimer](docs/Roadmpa.md#disclaimer)
- [Geometry](docs/Roadmap.md#geometry)
- [Solvers](docs/Roadmap.md#solvers)
- [Aerodynamics](docs/Roadmap.md#aerodynamics)
  - [Losses](docs/Roadmap.md#losses)
  - [Flows](docs/Roadmap.md#flows)
  - [Forces](docs/Roadmap.md#forces)
- [Other](docs/Roadmap.md#geometry)

---

## Disclaimer

As this code is not directly linked to my Ph.D. thesis, its development is not
my highest priority.

_The following features are listed in not specific order. Their development will
depend on my available time and current needs._

## Geometry

- [ ] Varying airfoil
- [x] Linear, custom or "ideal" twist and taper
- [ ] Fully customizable blade
- [ ] Coaxial rotors

## Solvers

- [x] Linearised equations for hovering rotors ([Leishman](docs/Solver-Leishman.md))
- [x] Full system with induction factors ([IndFact](docs/Solver-IndFact.md))
- [x] Full system with one single equation ([Stahlut](docs/Solver-Stahlut.md))
- [ ] Full system with alternative velocity triangles
- [ ] Improve stability of ([Stahlut](docs/Solver-Stahlut.md))
- [ ] Proper validation of all the solvers

## Aerodynamics

### Losses

- [x] Tip losses (Prandtl)
- [x] Hub losses (Prandtl)
- [ ] Spinner effects

### Flows

- [ ] Compressibility effects (use directly polar data at correct M?)
- [ ] Non-axial flows (helicopters in forward flight, tiltrotors, quadcopters)

### Forces

- [ ] Accept polynomial forms for $`C_l`$ and $`C_d`$ instead of using polars
- [ ] Integrate XFoil directly to recover polars on-the-fly

## Other

- [ ] Include [XFoil2Mat](https://gitlab.uliege.be/thlamb/xfoil2mat) as a
  submodule to directly process the XFoil data.
