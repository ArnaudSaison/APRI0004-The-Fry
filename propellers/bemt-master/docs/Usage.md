- [Installation](docs/Usage.md#installation)
- [Run the code](docs/Usage.md#run-the-code)
- [Customize the configuration](docs/Usage.md#customize-the-configuration)

---

## Installation

1. Ensures that Matlab is installed and the few [prerequisites](docs/Prerequisites.md)
   are met.
2. Clone the repository (or simply download it as a *.zip* file).

## Run the code

The code comes with a template configuration and a few airfoil polars. It is
possible to run it directly without any extra configuration.

1. Navigate to the `BEMT/src/` folder
2. Run `main.m`

## Customize the configuration

1. Make a copy of the template configuration
   `src/configurations/config_template.m`
   and modify it to specific the configuration you wish to study[^1].
2. Add the data for your airfoils in `src/airfoils_data/`. See [Airfoil
   data](docs/Code-AirfoilData.md) for more information.
4. Modify `main.m` to indicate your configuration file (l.60) and execute it.

[^1]: Making a copy allows you to always recover a proper configuration file in
  case you delete important parameters or comments.
