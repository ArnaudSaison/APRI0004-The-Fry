- [Matlab](docs/Prerequisites.md#matlab)
  - [Aerospace Toolbox](docs/Prerequisites.md#aerospace-toolbox)
- [XFoil](docs/Prerequisites.md#xfoil)
  - [XFoil2Mat](docs/Prerequisites.md#xfoil2mat)

---

## Matlab (mandatory)

_Mandatory requirement. Usage is not guaranteed on other platforms (such as
Octave)._

This code uses some Matlab functions that are not available in other systems
(such as Octave). The code was mainly developed using Matlab 2018a, it is not
guaranteed to run on older versions.

### Aerospace Toolbox

_Heavily recommended, requires workaround if not installed._

The code uses the Matlab function `atmosisa`, which is part of the _Aerospace
Toolbox_. This function simply calculates the air parameters with the
International Standard Atmosphere standard based on the altitude. It is used to
retrieve the air density. If this function is not present on your Maltab
installation, you can download an [alternate open-source
version][atmosisa-foss].

## XFoil

_Optional for testing, mandatory for proper analysis._

This code requires the loading of airfoil polars to determine the $`C_l`$ and
$`C_d`$ of each blade element. The common way to generate such polars is to use
[XFoil](https://web.mit.edu/drela/Public/web/xfoil/). Currently, XFoil is not
called automatically by the BEMT code, so the user must generate manually a few
polars for their airfoils of interest.

As stated in the [roadmap](docs/Roadmap.md), a future release may add direct
integration of XFoil, so the user does not have to manually generate the polars.

While XFoil is listed in the prerequisites, a few airfoil polars are already
bundeled in the present repository for testing and example purpose. So it is
perfectly possible to run this code without having XFoil installed.

### XFoil2Mat

_Optional for testing, mandatory for proper analysis._

The BEMT code does not process directly the XFoil data as it is. These raw XFoil
output need to be parsed and formatted in a clear structure first. To do that, a
separate script is provided in the form of
[XFoil2Mat](https://gitlab.uliege.be/thlamb/xfoil2mat).

This script will be integrated in the BEMT project in the future, in order to
streamline the processing of raw data for this code.

[atmosisa-foss]: https://github.com/lentzi90/octave-atmosisa
