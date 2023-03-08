from scipy.integrate import solve_ivp
import numpy as np

def trajectory(x0, func, T=1000, dt=0.01):
    
    ts = np.arange(0, T*dt, dt)
    sol = solve_ivp(lambda t, y: func(y), [0, T*dt], x0, t_eval=ts)
    
    return sol.y

def lorentz(x, sigma, rho, beta):
    x, y, z = x
    return np.array([sigma * (y - x), x * (rho - z) - y, x * y - beta * z])

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x0 = np.random.normal(size=3)
traj = trajectory(np.random.normal(size=3), lambda x: lorentz(x, 10, 28, 8.0/3), T=100000, dt=0.01)

print(traj.shape)