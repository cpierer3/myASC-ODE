# Implementing and Comparing ODE Methods

We implemented and compared three numerical methods for solving ordinary differential equations (ODEs).

## Numerical Methods

### Explicit Euler

$$
x_{n+1} = x_n + h\, f(t_n, x_n)
$$

### Implicit Euler

$$
x_{n+1} = x_n + h\, f(t_{n+1}, x_{n+1})
$$

### Improved Euler (Heun’s Method)

$$
x_{n+1} = x_n + \frac{h}{2}
\left(
f(t_n, x_n)
+ f\bigl(t_{n+1}, x_n + h\, f(t_n, x_n)\bigr)
\right)
$$

---

## Test Problem

To investigate the long-time stability of these methods, we solve the second-order ODE

$$
\ddot{y} = -y
$$

which describes a **simple harmonic oscillator**.

Rewriting it as a first-order system yields

$$
\dot{y} = v
$$

$$
\dot{v} = -y
$$

The exact solution conserves the total energy

$$
E = \frac{1}{2}\left(v^2 + y^2\right)
$$

making this problem well suited for studying numerical stability and energy conservation.

---

## Numerical Results

The following plots show the behavior of the numerical solutions for larger time intervals:

- **Phase portrait:** velocity vs. position  
- **Position vs. time**  
- **Velocity vs. time**

### Phase Portrait
![Phase portrait](../../demos/phase_plot.png)

### Position vs. Time
![Position vs Time](../../demos/position_vs_time.png)

### Velocity vs. Time
![Velocity vs Time](../../demos/velocity_vs_time.png)

---

## Discussion

The plots indicate that:

- The **explicit Euler method** adds energy to the system, leading to an unstable spiral in phase space.
- The **implicit Euler method** removes energy, causing the solution to decay over time.
- The **improved Euler method** conserves energy much better and closely reproduces the expected periodic motion.
