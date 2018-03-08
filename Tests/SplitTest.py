import numpy as np

a = np.array([10, 20, 0, 50, 60, 0, 70, 80, 0, 100, 20])

b = np.where(a == 0)[0]
print(b)

c = np.split(a, b)

print(c)