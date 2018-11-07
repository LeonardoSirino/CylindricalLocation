import numpy as np
from numba import jit
import time

values = np.linspace(-1, 1, num=2000000)


@jit(nopython=True)
def SomeSlowStuff(values):
    for value in values:
        x = np.arccos(value)
        return x


t0 = time.time()
x = SomeSlowStuff(values)
print(x)
t1 = time.time()
dt = t1 - t0

print("Primeira vez: " + str(dt))

t0 = time.time()
x = SomeSlowStuff(values)
print(x)
t1 = time.time()
dt = t1 - t0

print("Segunda vez: " + str(dt))

t0 = time.time()
for value in values:
    x = np.arccos(value)

t1 = time.time()
dt = t1 - t0

print("Sem otimizar: " + str(dt))
