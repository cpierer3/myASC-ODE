# myASC-ODE Simulation Framework

Welcome to the myASC-ODE simulation framework documentation. This book contains interactive simulations and analysis of mechanical systems using C++ and Python.

## Features

- **Kreisel Systems**: Simulating a spinning Kreisel
- **Crane Dynamics**: Realistic crane structure modeling with vibration analysis

## Interactive Simulations

This book contains fully interactive Jupyter notebooks with:
- Real-time 3D visualizations
- Dynamic parameter adjustments
- Comprehensive analysis tools

Only the user-friendly interface is implemented in Python; all performance-critical computations are handled by the underlying C++ framework.


Explore the interactive notebooks to see the simulations in action! We will first briefly mention some of the theory behind the implemented methods.

### Solving a Mass Spring System with a generalized $\alpha$ Method

In an effort to avoid the instabilities that the Newmark method leads to for non linear ODEs, we can introduce the generalized $\alpha$ method. We recall the Newmark scheme for a second order ODE

$$
M \ddot{x} = F(x)
$$

to be given as

$$
\begin{aligned}
x_{n+1} &= x_n + \tau v_n + \tau^2 \left( \left( \tfrac{1}{2} - \beta \right) a_n + \beta a_{n+1} \right), \\
v_{n+1} &= v_n + \tau \left( (1-\gamma) a_n + \gamma a_{n+1} \right)
\end{aligned}
$$

where $v_n = \dot{x}_n$ and $a_n = M^{-1} F(x_n)$.

For a damping parameter $\rho^\infty$, we can introduce the new variables

$$
\begin{aligned}
\alpha_m &= \frac{2 \rho^\infty - 1}{\rho^\infty + 1}, \\
\alpha_f &= \frac{\rho^\infty}{\rho^\infty + 1}
\end{aligned}
$$

and set

$$
\begin{aligned}
\beta &= \frac{1}{4} (1 - \alpha_m + \alpha_f)^2, \\
\gamma &= \frac{1}{2} - \alpha_m + \alpha_f
\end{aligned}
$$

to get to the generalized $\alpha$ method formulation

$$
\begin{aligned}
x_{n+1} &= x_n + \tau v_n + \tau^2 \left( \left( \tfrac{1}{2} - \beta \right) a_n + \beta a_{n+1} \right), \\
v_{n+1} &= v_n + \tau \left( (1-\gamma) a_n + \gamma a_{n+1} \right), \\
x_{n+1-\alpha_f} &= (1-\alpha_f) x_{n+1} + \alpha_f x_n, \\
a_{n+1-\alpha_m} &= (1-\alpha_m) a_{n+1} + \alpha_m a_n
\end{aligned}
$$

#### Systems with constraints

We want to implement a way to model joints between two masses. For this, we introduce systems with constraints. We define the Lagrange function for a constrained system as

$$
L(x, \lambda) = -U(x) + \langle \lambda, g(x) \rangle
$$

with

$$
U(x) = m g x_z, \qquad
g(x) := \lVert x - x_0 \rVert^2 - l^2
$$

being the length constraint.

This leads to the second order system of ODEs

$$
\begin{aligned}
m_i \ddot{x}_i &= \frac{\partial}{\partial x_i} L(x, \lambda), \\
0 &= \nabla_\lambda L(x, \lambda)
\end{aligned}
$$

This system can be solved using the generalized $\alpha$ method.
