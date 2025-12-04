import numpy as np
import matplotlib.pyplot as plt


import sys

if sys.platform.startswith("win"):
    path_jack     = r"C:\Users\lukas\Documents\Scicomp2\myASC-ODE\build\output_test_ode.txt"
    path_implicit = r"C:\Users\lukas\Documents\Scicomp2\myASC-ODE\build\output_test_ode_implicit.txt"
    path_explicit = r"C:\Users\lukas\Documents\Scicomp2\myASC-ODE\build\output_test_ode_explicit.txt"
    path_improved = r"C:\Users\lukas\Documents\Scicomp2\myASC-ODE\build\output_test_ode_improved.txt"
else:    
    path_jack     = r"output_test_ode.txt"
    path_implicit = r"output_test_ode_implicit.txt"
    path_explicit = r"output_test_ode_explicit.txt"
    path_improved = r"output_test_ode_improved.txt"

data_jack = None
data_implicit = None
data_explicit = None
data_improved = None

try:
    data_jack = np.loadtxt(path_jack, usecols=(0, 1, 2))
except Exception as e:
    print("Could not load jack:", e)

try:
    data_implicit = np.loadtxt(path_implicit, usecols=(0, 1, 2))
except Exception as e:
    print("Could not load implicit:", e)

try:
    data_explicit = np.loadtxt(path_explicit, usecols=(0, 1, 2))
except Exception as e:
    print("Could not load explicit:", e)

try:
    data_improved = np.loadtxt(path_improved, usecols=(0, 1, 2))
except Exception as e:
    print("Could not load improved:", e)

# Position vs Time
plt.figure()
if data_jack is not None:
    plt.plot(data_jack[:,0], data_jack[:,1], label="jack")
if data_implicit is not None:
    plt.plot(data_implicit[:,0], data_implicit[:,1], label="implicit")
if data_explicit is not None:
    plt.plot(data_explicit[:,0], data_explicit[:,1], label="explicit")
if data_improved is not None:
    plt.plot(data_improved[:,0], data_improved[:,1], label="improved")
plt.xlabel("time")
plt.ylabel("position")
plt.title("Position vs Time")
plt.grid(True)
plt.legend()

# Velocity vs Time
plt.figure()
if data_jack is not None:
    plt.plot(data_jack[:,0], data_jack[:,2], label="jack")
if data_implicit is not None:
    plt.plot(data_implicit[:,0], data_implicit[:,2], label="implicit")
if data_explicit is not None:
    plt.plot(data_explicit[:,0], data_explicit[:,2], label="explicit")
if data_improved is not None:
    plt.plot(data_improved[:,0], data_improved[:,2], label="improved")
plt.xlabel("time")
plt.ylabel("velocity")
plt.title("Velocity vs Time")
plt.grid(True)
plt.legend()

# Phase plot: velocity vs position
plt.figure()
if data_jack is not None:
    plt.plot(data_jack[:,1], data_jack[:,2], label="jack")
if data_implicit is not None:
    plt.plot(data_implicit[:,1], data_implicit[:,2], label="implicit")
if data_explicit is not None:
    plt.plot(data_explicit[:,1], data_explicit[:,2], label="explicit")
if data_improved is not None:
    plt.plot(data_improved[:,1], data_improved[:,2], label="improved")
plt.xlabel("position")
plt.ylabel("velocity")
plt.title("Phase Plot (v vs x)")
plt.grid(True)
plt.legend()

plt.show()