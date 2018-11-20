import numpy as np

a = np.array([1, 2, 3])

b = a - np.min(a)
b = np.exp(-b)

print(b)