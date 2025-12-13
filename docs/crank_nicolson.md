# Crank-Nicolson Method

The Crank-Nicolson method is a second-order implicit numerical method that combines the best features of explicit and implicit approaches. It's particularly renowned for its excellent stability properties and is widely used in scientific computing.

## Mathematical Formulation

Given an initial value problem:
$$\frac{dy}{dt} = f(t, y), \quad y(t_0) = y_0$$

The Crank-Nicolson method uses the trapezoidal rule:
$$y_{n+1} = y_n + \frac{h}{2}[f(t_n, y_n) + f(t_{n+1}, y_{n+1})]$$

This averages the derivatives at both the current and next time steps.

## Algorithm

1. Start with initial conditions: $t_0, y_0$
2. Choose a step size $h$
3. For each step $n = 0, 1, 2, ...$:
   - Evaluate $f_n = f(t_n, y_n)$
   - Solve the implicit equation: $y_{n+1} - y_n - \frac{h}{2}[f_n + f(t_{n+1}, y_{n+1})] = 0$
   - Update: $t_{n+1} = t_n + h$

