import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import time

T = 1
A = 1
N = 50
w = 2 * np.pi / T

@jit(nopython=True, parallel=True)
def x_t(t, N):
    x = 0
    for n in range(1, N):
        bn = 2 * A / (np.pi * n)
        x += bn * np.sin(n * w * t)
    
    return x

size = 10000
t = np.linspace(0, 1, num=size)
legend = []
t0 = time.time()
for N in np.linspace(10, 50, num=5):
    x_vec = np.zeros(size)
    N = int(N)
    k = 0
    for temp in t:
        x_vec[k] = x_t(temp, N)
        k += 1

    plt.plot(t, x_vec)
    legend.append(str(N))

t1 = time.time()
dt = t1 - t0
print(dt)
plt.legend(legend)
plt.show()