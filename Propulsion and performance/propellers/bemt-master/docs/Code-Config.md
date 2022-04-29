- [Parameters](docs/Code-Config.md#parameters)
  - [Simulation](docs/Code-Config.md#simulation)
  - [Free stream](docs/Code-Config.md#free-stream)
- [Blades](docs/Code-Config.md#blades)
- [Elements](docs/Code-Config.md#elements)

---

The configuration file is divided in different segments, each one corresponding
to a different Matlab structure.

A template configuration file can be found in
[src/configuration/config_template](https://gitlab.uliege.be/thlamb/bemt/-/blob/master/src/configurations/config_template.m).

## Parameters

### Simulation

The following parameter controls the simulation in itself. It focuses on the
technical details such as display, saves, plots, etc. It also controls the
simulation behaviors and some models (losses, convergence criterion, etc.)

**Structure name**: `Param.Simul`

| **Field**     | **Use**                     | **Type**   | **Example**               |
| SAVE_RESULTS  | Toggle autosave of results  | logical    | `true`, `false`, `0`, `1` |
| SAVE_PATH     | Directory for saved results | string     | `results/`                |
| SAVE_FILENAME | Filename of saved results   | string     | `myrotor`                 |
| SHOW_GRAPHS   | Toggle display of plots     | logical    | `true`, `false`, `0`, `1` |
| SHOW_3DVIEW   | Toggle display of 3D view   | logical    | `true`, `false`, `0`, `1` |
| PRINT_CONSOLE | Toggle console output       | logical    | `true`, `false`, `0`, `1` |
|  |  |  |  |
| SYSTEM_SOLVER | Select solver               | string | `Leishman`,`IndFact`,`Stahlut` |
| LOSSES        | Include Prandtl losses      | string | `none`, `tip`, `root`, `all`  |
| TWIST         | Twist model                 | string | `linear`, `ideal`, `custom`   |
| TAPER         | Taper model                 | string | `linear`, `ideal`, `custom`   |
| TRIM_IDEAL    | Trim ideal twist and taper to realistic values | logical | `true`, `false`, `0`, `1` |
| CONV_CRIT     | Convergence criterion       | scalar | - |

### Free stream

These parameter set the values related of the free stream.

**Structure name**: `Param.Air`

| **Field**     | **Use**                      | **Type** | **Unit** |
| ------------- | ---------------------------- | -------- | -------- |
| ALTITUDE      | Flight altitude              | scalar   | m        |
| AXIAL_VELOC   | Freestream axial velocity    | scalar   | m/s      |
| DYN_VISCOSITY | Dynamic viscosity of the air | scalar   | Pa.s     |

## Blades

General geometric parameters of a blade.

**Structure name**: `Blades`

| **Field**     | **Use**                                    | **Type** | **Unit** |
| ------------- | ------------------------------------------ | ------   | ---      |
| PROFILE_FILE  | `.dat` file with the profile's coordinates | string   |          |
| COEFF_FILE    | `.mat` file with the airfoil polars        | string   |          |
| OMEGArpm      | Rotational speed of the propeller          | scalar   | RPM      |
| nELEM         | Number of blade elements                   | scalar   | -        |
| nBLADES       | Number of blades in the rotor              | scalar   | -        |
| RADIUS        | Rotor radius (from hub center to tip)      | scalar   | m        |
| ROOT_CHORD    | Blade root chord                           | scalar   | m        |
| TIP_CHORD     | Blade tip chord                            | scalar   | m        |
| ROOT_CUTOUT   | Blade root cutout                          | scalar   | m        |
| THETA_TIPdeg  | Stagger angle at the tip                   | scalar   | deg      |
| COLL_PITCHdeg | Rotor collective pitch angle               | scalar   | deg      |

## Elements

_Optional. Used to specify custom, non-linear, twist or taper._

**Structure name**: `Elem`

| **Field**    | **Use**                       | **Type**           | **Unit** |
| ------------ | ----------------------------  | ------------------ | ---- |
| CUSTOM_TWIST | Stagger angle of each element | vector `[1,nELEM]` | deg  |
| CUSTOM_CHORD | Chord of each element         | vector `[1,nELEM]` | m    |
