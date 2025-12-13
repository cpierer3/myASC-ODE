# Implicit Euler Method

The Implicit Euler method is a first-order implicit numerical method for solving ordinary differential equations. Unlike the explicit version, it evaluates the derivative at the new time step, providing better stability properties.

## Mathematical Formulation

Given an initial value problem:
$$\frac{dy}{dt} = f(t, y), \quad y(t_0) = y_0$$

The Implicit Euler method uses:
$$y_{n+1} = y_n + h \cdot f(t_{n+1}, y_{n+1})$$

This creates an implicit equation that must be solved for $y_{n+1}$ at each step.

