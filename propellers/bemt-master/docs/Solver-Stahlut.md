- [Solver](docs/Solver-Stahlut.md#solver)
- [Hypotheses](docs/Solver-Stahlut.md#hypotheses)
- [Limitations](docs/Solver-Stahlut.md#limitations)

---

## Solver

The `Stahlut` solver relies on the recent work of Stahlut and Leishman, where
they were able to express the whole system as one single transcendental
equation.

Thanks to that formulation, this solver is more precise that the
[Leishman](docs/Solver-Leishman.md), which relies on numerous extra assumptions.
It is also often more stable than the [IndFact](docs/Solver-IndFact.md) solver
and does not come with the limitation related to the airspeed.

| **Appropriate for**     | Everything  |
| **Not appropriate for** |             |

## Hypotheses

- This solver does not require any additional hypotheses (besides the common
  hypotheses for the blade element theory).

## Limitations

- More time consuming.
- Resolution with `fzero` is not always guaranteed

[^1]: Stahlhut, C. and Leishman, J.G. *Aerodyanamic Design Optimization of
  Proprotors for Convertible-Rotor Concepts*, American Helicopter Society 68th
  Annual Forum, 1-3 May 2012.
