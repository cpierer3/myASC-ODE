# myASC-ODE Simulation Framework

Welcome to the myASC-ODE simulation framework documentation. This book contains interactive simulations and analysis of mechanical systems using C++ and Python.

## Features

- **Mass-Spring Systems**: Comprehensive 3D mass-spring system simulations
- **Crane Dynamics**: Realistic crane structure modeling with vibration analysis
- **Kreisel Systems**: Gyroscopic and rotational dynamics simulations

## Interactive Simulations

This book contains fully interactive Jupyter notebooks with:
- Real-time 3D visualizations
- Dynamic parameter adjustments
- Comprehensive analysis tools
- C++ and Python implementations  


Explore the interactive notebooks to see the simulations in action! We will first briefly mention some of the theory behind the implemented methods.

### Solving a Mass Spring System with a generalized $\alpha$ Method

In an effort to avoid the instabilities that the newmark method leads to for non linear ODEs, we can Introduce the generalized $\alpha$ method. We recall the Newmark Scheme for a second order ODE $ M\ddot x = Fx $ to be given as
$$
x_{n+1} & = & x_n + \tau v_n  + \tau^2 ( ( \frac{1}{2} - \beta) a_n + \beta a_{n+1} ) \\
v_{n+1} & = & v_n + \tau ( (1-\gamma) a_n + \gamma a_{n+1} )
$$
where $v_n$ is $\dot x_n$ and $a_n = M^{-1}F(x_n)$. For a damping parameter $ \rho^\infty $ We can introduce the new variables 
\begin{align*}
\alpha_m & = & \frac{2 \rho^\infty - 1}{\rho^\infty + 1} \\
\alpha_f & = & \frac{\rho^\infty}{\rho^\infty + 1} \\
\end{align*}
and set 
\begin{align*}
\beta & = & \frac{1}{4} (1 - \alpha_m + \alpha_f)^2 \\
\gamma & = & \frac{1}{2} - \alpha_m + \alpha_f
\end{align*}
to get to the generalized $\alpha$ method formulation
\begin{align*}
x_{n+1} & = & x_n + \tau v_n  + \tau^2 ( ( \frac{1}{2} - \beta) a_n + \beta a_{n+1} ) \\
v_{n+1} & = & v_n + \tau ( (1-\gamma) a_n + \gamma a_{n+1} ) \\
x_{n+1-\alpha_f} & = & (1-\alpha_f) x_{n+1} + \alpha_f x_n \\
a_{n+1-\alpha_m} & = & (1-\alpha_m) a_{n+1} + \alpha_m a_n
\end{align*}

#### Systems with constraints

We want to implement a way to model joints between two masses, for this we introduce systems with constraints. We can define the Lagrange function for a system with constraints as
$$
L(x, \lambda) = -U(x) + \left< \lambda , g(x) \right>,
$$
with $U(x) = m g x_z$ and $g(x) := | x - x_0 |^2 - l^2$ being the length constraint, and derive the second Order System of ODEs
$$
m_i \ddot x_i  =  \frac{\partial}{\partial x_i} L(x, \lambda) \\
 0  =  \nabla_\lambda L(x, \lambda)
$$
This we can solve using the generalized $\alpha$ method. 
