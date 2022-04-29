# BEMT

**BEMT** is an implementation of the Blade Element Momentum Theory in Matlab,
used for teaching purposes at the [University of
Liege](https://www.am.uliege.be/) (Belgium).  It allows students to quickly test
different rotors/propellers configurations and see how they perform in terms of
Thrust, Power and Torque.

This code does not accommodate varying airfoil, but it can take into account any
twist or taper. It can also calculate the ideal twist or taper using analytical
formulas, based on some hypotheses. The tip losses can be modeled using Prandlt
tip loss function.

The current implementation is only suited for single rotors/props in isolation.
Three different solvers are implemented: one follows the methodology commonly
used for propeller design (itertive scheme with inflow factors _a_ and _b_), one
is a linearized solution as described by Leishman (more suited for helicopters)
and the third one uses a more complete formulation as described by Stalhut. See
the [Solvers
documentation](https://gitlab.uliege.be/thlamb/bemt/wikis/BEMT/Solvers) for more
details.

## Installation and utilization

1. Clone the repository (or simply download it as a *.zip* file).
2. Create a new config file based on `config_template.m` and modify it with your
   rotor/propeller configuration.
3. Add the data for your airfoils in `src/airfoils_data`. See [Airfoil
   data](https://gitlab.uliege.be/thlamb/bemt/wikis/BEMT/Airfoil%20data)
   documentation.
4. Modify `main.m` to indicate your configuration file (l.60) and execute it.

## Documentation

See [the wiki][code_wiki] for a detailed explanation of the code.

The slides of the presentation are in a [separate repository][slides_repo]. They
give a theoretical explanation of the BEMT and a general overview of the solving
process.  You can use the following links to directly download the latest
version

- For helicopters: [Slides][slides_heli] - [Handout][handout_heli]
- For propellers: [Slides][slides_props] - [Handout][handout_props]

[code_wiki]: https://gitlab.uliege.be/thlamb/bemt/wikis/home
[slides_repo]:https://gitlab.uliege.be/thlamb/slides-and-presentations
[slides_heli]:https://thlamb.gitlabpages.uliege.be/slides-and-presentations/bemt_helico_slides.pdf
[handout_heli]:https://thlamb.gitlabpages.uliege.be/slides-and-presentations/bemt_helico_handout.pdf
[slides_props]:https://thlamb.gitlabpages.uliege.be/slides-and-presentations/bemt_prop_slides.pdf
[handout_props]:https://thlamb.gitlabpages.uliege.be/slides-and-presentations/bemt_prop_handout.pdf
