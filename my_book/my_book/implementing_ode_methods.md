## Implementing ODE Methods
We have implemented and compared three different mehtods to solve ODEs:
Explicit Euler: $$x_{n+1} = x_n + h\, f(t_n, x_n)$$

Implicit Euler: $$x_{n+1} = x_n + h\, f(t_{n+1}, x_{n+1})$$

Improved Euler: $$x_{n+1} = x_n + \tfrac{h}{2}\,\big[f(t_n, x_n) + f\big(t_{n+1}, x_n + h\, f(t_n, x_n)\big)\big]$$

In the following images we want to see whether the methods are stable for larger times. We are solving the second order ODE defined by $\ddot y = -y$. On the first and Second plots we the the positions plotted against time, on the third plot we see the velocity against time. The plots suggest that that the explicit euler add energy to the system, while the implicit euler removes energy from the system. The improved euler method seems to conserve energy quite well. 

![Phasenportrait](../../demos/phase_plot.png)
![Position vs Time](../../demos/position_vs_time.png)
![Velocity vs Time](../../demos/velocity_vs_time.png)
