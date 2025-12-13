# Improved Euler Method (Heun's Method)

The Improved Euler method, also known as Heun's method, is a second-order Runge-Kutta method that provides significantly better accuracy than the standard Euler methods while maintaining reasonable computational cost.

## Mathematical Formulation

Given an initial value problem:
$$\frac{dy}{dt} = f(t, y), \quad y(t_0) = y_0$$

The Improved Euler method uses a two-stage process:
$$k_1 = f(t_n, y_n)$$
$$k_2 = f(t_n + h, y_n + h \cdot k_1)$$
$$y_{n+1} = y_n + \frac{h}{2}(k_1 + k_2)$$
