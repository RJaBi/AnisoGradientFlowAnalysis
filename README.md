# AnisoGradientFlowAnalysis
A repository to do the analysis of the Gradient Flow output files from the ancillary files of arXiv: 1203.4469 for the purposes of gauge anisotropy determination and scale setting (lattice spacing determination).

Specifically the code takes a `toml` input file detailing the various parameters and data files, determines the anisotropy using the $R_E$ method as in arXiv: 1205.0781 and then determines the (spatial) lattice spacing using the $W_0$ scale approach. At all stages, [b-splines](https://github.com/jacobwilliams/bspline-fortran) are used to interpolate between the flow anisotropy values.



This repository uses the [fortran-lang/fpm](https://github.com/fortran-lang/fpm) to build. See the instructions below on how to get this. This package manager automatically gets the quired dependencies from GitHub. If there is no internet connection or git available you will need acquire them manually and then modify the `fpm.toml` file to i.e.
```
[dependencies]
stdlib = { path = 'myrelativepathto/stdlib'}
```


## How to install the FPM

See the instructions [here](https://fpm.fortran-lang.org/install/index.html)
If you are on a Linux distro and have any Fortran compiler installed, do:
```
git clone https://github.com/fortran-lang/fpm
cd fpm
./install.sh
```

This will put the fpm in your `$HOME/.local/bin`

## How to change the FPM config

To build, simply: `fpm build`


To run: `fpm run -- FULL_PATH_TO_TOML_FILE`



## pre-commit hooks

The repo also comes with a pre-commit that will ensure a formatting for your Fortran files. You can install pre-commits by using:

```
python -m pip install pre-commit
pre-commit install
```
