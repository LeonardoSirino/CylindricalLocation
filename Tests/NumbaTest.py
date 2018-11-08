from numba import jit
import time
import numpy as np

values = np.linspace(0, 1, num=1000000)


@jit(nopython=True)
def SomeSlowStuff(values):
    for value in values:
        x = np.arccos(value)

    return x


def SomeSlowerStuff(values):
    for value in values:
        x = np.arccos(value)

    return x


t0 = time.time()
x = SomeSlowStuff(values)
t1 = time.time()
dt = t1 - t0

print("Primeira vez: " + str(dt))

t0 = time.time()
x = SomeSlowStuff(values)
t1 = time.time()
dt = t1 - t0

print("Segunda vez: " + str(dt))

t0 = time.time()
x = SomeSlowerStuff(values)
t1 = time.time()
dt = t1 - t0

print("Sem otimizar: " + str(dt))
