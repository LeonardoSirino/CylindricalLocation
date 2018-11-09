import numpy as np
from numba import jit

@jit(nopython=True, parallel=True)
def wallDist(x1, y1, x2, y2, d):
    dist1 = np.sqrt((x1 - x2) ** 2 + (y1 - y2)**2)
    # Clone à direita
    dist2 = np.sqrt((x1 - x2 + d * np.pi) ** 2 + (y1 - y2)**2)
    # Clone à esquerda
    dist3 = np.sqrt((x1 - x2 - d * np.pi) ** 2 + (y1 -y2)**2)

    if dist1 < dist2:
        dist = dist1
    else:
        dist = dist2

    if dist3 < dist:
        dist = dist3
    
    return dist