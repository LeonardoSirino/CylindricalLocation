import matplotlib.pyplot as plt
import math as m

C = 2492.0
d = C / m.pi
h = 2700.0
f = 0.5
sp = 510

x_vessel = [0, C, C, 0, 0, 0, C, C, 0, 0, C, C, 0, 0]
y_vessel = [0, 0, h, h, 0, -sp, -sp, 0, 0, h, h, h + sp, h + sp, h]

plt.plot(x_vessel, y_vessel)
plt.show()