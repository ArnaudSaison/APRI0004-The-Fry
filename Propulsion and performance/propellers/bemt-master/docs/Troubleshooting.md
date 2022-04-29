- [Unknown function atmosisa](docs/Troubleshooting.md#unknown-function-atmosisa)
- [fzero crashes](docs/Troubleshooting.md#fzero-crashes)
- [Error when using IndFact solver](docs/Troubleshooting.md#error-when-using-indfact-solver)
- [LinearRange:BoundNotFound](docs/Troubleshooting.md#linearrangeboundnotfound)
- [Other issues](docs/Troubleshooting.md#other-issues)

## Unknown function atmosisa

As explained in the [Prerequisites](docs/Prerequisites.md#Aerospace-toolbox),
this code uses Matlab's function `atmosisa` (from the Aerospace Toolbox) in
order to get the ISA values for a given altitude.

If this function is unavailable to you, please download a free and open source
variant it from [this Github
repository](https://github.com/lentzi90/octave-atmosisa) and place it somewhere
in Matlab's Path (either in the root of this project, or in `Documents/MATLAB`)

## fzero crashes

Matlab's function `fzero` is used in the [Stahlut](docs/Solver-Stahlut.md)
solver in order to find the roots of the transcendental equation.

This function relies on the bracketed bisection method, which is initiated with
bouds for the angle of incidence at 0° and 90°. In some rare circumstances, the
evaluation of the function in these two bounds may return values that have the
same sign. If this is the case, *fzero* will crash as it can not be sure that a
root lies between these bounds. Lowering the upper bound should normally solve
this issue. If this does not work, please submit an issue report detailing the
problem and Matlab's error message.

_Improvements related to the stability of this solver are already part of
the [roadmap](docs/Roadmap.md)_.

## Error when using `IndFact` solver

The [IndFact](docs/Solver-IndFact.md) solver solves the full system of equations
with an iterative process. As the system of equation is nonlinear, the
convergence of the iterative process is not guaranteed. To mitigate that, a
relaxation factor is introduced in the solver, so the loops should move slowly
but steadily towards to correct solution.

If that is not the case, some user adjustments may be required to make the
solver work with the specific conditions and geometry.

The most frequent issues are:

- Tip loss factor is too small at the last blade element, which leads to a
  singular expression for the inflow factors. A solution would be to bound the
  minimum tip loss to half the one of the previous element for instance.
- Solver does not converge or other errors. In order to facilitate the
  convergence, a relaxation factor is defined in this solver function. Lowering
  it a bit more might help with convergence, at the expense of computation time.

Currently these tweaks need to be hard-coded by the user facing the issue.

In any case, this solver solves the exact same equations base equations than the
[Stahlut](docs/Solver-Stahlut.md) solver. It is therefore advise to use it
instead of this one if you are facing some issues.

_Improvements related to the stability of this solver are already part of
the [roadmap](docs/Roadmap.md)_.

## LinearRange:BoundNotFound

In order to solve the BEMT, some specific values need to be extracted from the
airfoil polars (Lift curve specifically). It includes the zero-lift angle, the
stall angle and the lift-curve slope. In order to find the lift-curve slop. To
do that, a function needs to determine the linear range of the $`C_l-\alpha`$
curve by using an approximation of the first and second derivative. In some
rare cases, the polar does not include enough points or includes particular
conditions (deep stall) that lead to the impossibility of finding suitable
bounds for the linear range.

This error can be bypassed by changing `forceLinResults` in the function
`calcPolarChar`, although this is not advised. In this case, the linear range
bound will be replaced by the maximal/minimal angle of attack available. Keep in
mind that this will lead to incorrect results for the lift curve slope and
therefore for the overall outcome of the calculation.

If the Polar does not include any particular element (no double peaks, etc), a
solution might be to tweak the calculation parameters in the "Constants" section
of `calcPolarChar`.  If the polar is not regular, the proper way to troubleshoot
this would be to regenerate the polars using more suitable angles of attack in
order to limit the data to common operating conditions.

## Other issues

If you encounter any other issue with this code, please check first [the issue
tracker](https://gitlab.uliege.be/thlamb/bemt/issues) and fill a new issue
report if applicable. You can also contact me directly at tlambert@uliege.be.
