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

## Available Time-Stepping Methods

This library implements several numerical methods for solving ordinary differential equations. Each method has different properties regarding accuracy, stability, and computational cost. The following time-stepping methods are available:

1. **Explicit Euler Method** - A simple first-order explicit method
2. **Implicit Euler Method** - A first-order implicit method with better stability
3. **Improved Euler Method** - A second-order explicit method (also known as Heun's method)
4. **Crank-Nicolson Method** - A second-order implicit method with excellent stability properties

Each method is suitable for different types of problems depending on the characteristics of the differential equation system you want to solve.

## Method Details

### Explicit Euler Method

The Explicit Euler method is probably the simplest numerical method for solving ODEs. It follows:

$$y_{n+1} = y_n + h \cdot f(t_n, y_n)$$

**Characteristics:**
- **Order of accuracy:** First-order (O(h))
- **Stability:** Conditionally stable


### Implicit Euler Method

The Implicit Euler method improves stability by using an implicit formulation:

$$y_{n+1} = y_n + h \cdot f(t_{n+1}, y_{n+1})$$

**Characteristics:**
- **Order of accuracy:** First-order (O(h))
- **Stability:** Unconditionally stable for linear problems
- **Computational cost:** Higher (requires solving nonlinear equations)
- **Best for:** Stiff differential equations


### Improved Euler Method (Heun's Method)

The Improved Euler method is a second-order Runge-Kutta method that provides better accuracy:

$$k_1 = f(t_n, y_n)$$
$$k_2 = f(t_n + h, y_n + h \cdot k_1)$$
$$y_{n+1} = y_n + \frac{h}{2}(k_1 + k_2)$$

**Characteristics:**
- **Order of accuracy:** Second-order (O(h²))
- **Stability:** Conditionally stable
- **Computational cost:** cpmared to standard euler, it requires two function evaluations

### Crank-Nicolson Method

The Crank-Nicolson method combines implicit and explicit approaches for optimal stability and accuracy:

$$y_{n+1} = y_n + \frac{h}{2}[f(t_n, y_n) + f(t_{n+1}, y_{n+1})]$$

**Characteristics:**
- **Order of accuracy:** Second-order (O(h²))
- **Stability:** Unconditionally stable
- **Computational cost:** Higher (requires solving nonlinear equations)




   
