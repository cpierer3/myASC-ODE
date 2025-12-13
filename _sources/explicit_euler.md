# Explicit Euler Method

The Explicit Euler method is the simplest numerical method for solving ordinary differential equations (ODEs). It's a first-order method that provides a basic introduction to numerical ODE solving.

## Mathematical Formulation

Given an initial value problem:
$$\frac{dy}{dt} = f(t, y), \quad y(t_0) = y_0$$

The Explicit Euler method approximates the solution using:
$$y_{n+1} = y_n + h \cdot f(t_n, y_n)$$

## Algorithm

1. Start with initial conditions: $t_0, y_0$
2. Choose a step size $h$
3. For each step $n = 0, 1, 2, ...$:
   - Compute $f(t_n, y_n)$
   - Update: $y_{n+1} = y_n + h \cdot f(t_n, y_n)$
   - Update: $t_{n+1} = t_n + h$

## Properties

### Advantages
- **Simple implementation:** Only requires one function evaluation per step
- **Low computational cost:** Minimal memory and processing requirements
- **Easy to understand:** Straightforward mathematical concept

### Disadvantages
- **Low accuracy:** First-order method (error proportional to $h$)
- **Stability issues:** Can become unstable for large step sizes
- **Conditional stability:** Step size must be small enough to maintain stability

## Implementation in myASC-ODE

*[Add specific details about how to use the ExplicitEuler class in your library]*

```cpp
// Example usage (to be filled in with actual API)
ExplicitEuler solver;
// Implementation details here
```

## Examples



## When to Use

The Explicit Euler method is best suited for:
- Educational purposes and understanding basic numerical methods
- Simple problems with smooth solutions
- Cases where computational resources are extremely limited
- Initial prototyping before implementing more sophisticated methods

For production use, consider higher-order methods like Improved Euler or Runge-Kutta methods for better accuracy.
