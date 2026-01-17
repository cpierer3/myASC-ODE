# Welcome to myASC-ODE's documentation!

myASC-ODE is a C++ library for solving ordinary differential equations (ODEs).
The equation is defined by the right hand side function.
ASC-ODE provides various time-steppers which may be used for ODEs with right hand sides
given by a function object.

## Installation

Install myASC-ODE via git-clone:

```bash
git clone https://github.com/cpierer3/myASC-ODE.git
```

To configure and build some tests do:

```bash
cd my-ode-solver
mkdir build
cd build
cmake ..
make
```