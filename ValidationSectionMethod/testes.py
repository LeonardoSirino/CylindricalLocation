import math as m
import numpy as np
import matplotlib.pyplot as plt

b = 4
for a in np.linspace(0.1, 10, num=10):
    x_vec = np.linspace(-1, 1, num=500)
    y_vec = []
    for x in x_vec:
        y = a * np.log(0.5 * (x + 1)/(1-0.5*(x+1))) + b
        y_vec.append(y)

    plt.plot(x_vec, y_vec)

plt.show()