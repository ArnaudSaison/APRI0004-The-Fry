- [Solver](docs/Solver-IndFact.md#solver)
- [Hypotheses](docs/Solver-IndFact.md#hypotheses)
- [Limitations](docs/Solver-IndFact.md#limitations)

---

## Solver

The `IndFact` solver relies the most common BEMT description for propellers.
In that methodology, the velocities at the disk are expressed with the
deifnition of _induction factors_ (often called _a_ for the axial induction
factor and _b_ for the swirl factor):

```math
\left\lbrace\begin{array}{ll}V_{ax} &= V_\infty (1+a) \\ V_\theta &= \Omega y
(1-b)\end{array}\right.
```

Contrary to the [Leishman](docs/Solver-Leishman.md) implementation, this
formulation does not require extra hypotheses. It has nonetheless some
drawbacks. The most important one is that, due to the way the velocities are
defined, this solver can not be used with an airspeed of 0 (_i.e._ propeller at
start or hovering rotor).
This solver is also slower than [Leishman](docs/Solver-Leishman.md) and its
convergence is not always guaranteed.

| **Appropriate for**     | Everything but a few cases              |
| **Not appropriate for** | Hovering rotors and propellers at start |

## Hypotheses

- This solver does not require any additional hypotheses (besides the common
  hypotheses for the blade element theory).

## Limitations

- This solver can not treat the cases where the airspeed is 0 (hovering rotors).
- The convergence is not always guaranteed and may depend on various numerical
  parameters.

These two limitations can be alleviated by using the
[Stahlut](docs/Solvers-Stahlut.md) solver.
