- [Solver](docs/Solver-Leishman.md#solver)
- [Hypotheses](docs/Solver-Leishman.md#hypotheses)
- [Limitations](docs/Solver-Leishman.md#limitations)

---

## Solver

The `Leishman` solver relies on a linearization of the BEMT classical equations.
The methodology followed is described in full in Leishman's _Principle of
Helicopter Aerodynamics_[^1].

This methodology is often used by people interested in the analysis of hovering
rotors.

The linearization of equations introduces _de facto_ a few extra hypotheses and
limitations to this solver. Although, it comes with the benefit of being
extremely fast and with a guaranteed convergence.

| **Appropriate for**     | Helicopters and drones in hover                           |
| **Not appropriate for** | Propellers, Wind Turbines, large axial/forward velocities |

## Hypotheses

The linearization of equations lead to the following hypotheses:

- Small angle approximation
- Neglect $`dD`$ in front of $`dL`$
- Neglect swirl velocity

## Limitations

Due to the hypotheses mentioned above, the range of reliability is very limited.
Typically only hovering rotors or rotors in slow axial flow can be studied with
a decent validity.

The solve does not have any specific numerical limitations. It can be used with
any airspeed or condition (but then its reliability may not be guaranteed).

[^1]: Leishman, _Principles of Helicopter Aerodynamics_, Cambridge University
Press, 2006
