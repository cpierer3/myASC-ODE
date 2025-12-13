# Explicit Euler Method

The Explicit Euler method is the simplest numerical method for solving ordinary differential equations (ODEs). It's a first-order method that provides a basic introduction to numerical ODE solving.

## Mathematical Formulation

Given an initial value problem:
$$\frac{dy}{dt} = f(t, y), \quad y(t_0) = y_0$$

The Explicit Euler method approximates the solution using:
$$y_{n+1} = y_n + h \cdot f(t_n, y_n)$$
