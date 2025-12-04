import numpy as np
import matplotlib.pyplot as plt
import sys

if sys.platform.startswith("win"):
    path_crank = r"C:\Users\lukas\Documents\Scicomp2\myASC-ODE\build\output_test_ode_cranck.txt"
    path_improved = r"C:\Users\lukas\Documents\Scicomp2\myASC-ODE\build\output_test_ode_improved.txt"
else:
    path_crank = r"output_test_ode_crank.txt"
    path_improved = r"output_test_ode_improved.txt"

data_crank = None
data_improved = None

try:
    data_crank = np.loadtxt(path_crank, usecols=(0, 1, 2))
except Exception as e:
    print("Could not load Crank Nicolson data:", e)

try:
    data_improved = np.loadtxt(path_improved, usecols=(0, 1, 2))
except Exception as e:
    print("Could not load Improved Euler data:", e)

if data_crank is None and data_improved is None:
    raise SystemExit("No data loaded. Check file paths.")

plt.figure()
if data_crank is not None:
    plt.plot(data_crank[:,0], data_crank[:,1], label="Crank Nicolson")
if data_improved is not None:
    plt.plot(data_improved[:,0], data_improved[:,1], label="Improved Euler")
plt.xlabel("time")
plt.ylabel("position")
plt.title("Position vs Time")
plt.grid(True)
plt.legend()

# Velocity vs Time
plt.figure()
if data_crank is not None:
    plt.plot(data_crank[:,0], data_crank[:,2], label="Crank Nicolson")
if data_improved is not None:
    plt.plot(data_improved[:,0], data_improved[:,2], label="Improved Euler")
plt.xlabel("time")
plt.ylabel("velocity")
plt.title("Velocity vs Time")
plt.grid(True)
plt.legend()

# Phase plot: velocity vs position
plt.figure()
if data_crank is not None:
    plt.plot(data_crank[:,1], data_crank[:,2], label="Crank Nicolson")
if data_improved is not None:
    plt.plot(data_improved[:,1], data_improved[:,2], label="Improved Euler")
plt.xlabel("position")
plt.ylabel("velocity")
plt.title("Phase Plot (v vs x)")
plt.grid(True)
plt.legend()

plt.show()