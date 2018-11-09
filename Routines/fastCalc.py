import numpy as np
from numba import jit


@jit(nopython=True, parallel=True)
def wallDist(x1, y1, x2, y2, d):
    dist1 = np.sqrt((x1 - x2) ** 2 + (y1 - y2)**2)
    # Clone à direita
    dist2 = np.sqrt((x1 - x2 + d * np.pi) ** 2 + (y1 - y2)**2)
    # Clone à esquerda
    dist3 = np.sqrt((x1 - x2 - d * np.pi) ** 2 + (y1 - y2)**2)

    if dist1 < dist2:
        dist = dist1
    else:
        dist = dist2

    if dist3 < dist:
        dist = dist3

    return dist


@jit(nopython=True, parallel=True)
def sectionPos(s, pol):
    r = 1
    y = 0
    for coef in pol:
        y += coef * r
        r *= s

    return y

@jit(nopython=True, parallel=True)
def sectionArc(ri, rf, pol):
    x1 = 1
    x2 = 1
    s1 = 0
    s2 = 0
    for coef in pol:
        s1 += coef * x1
        x1 *= ri
        s2 += coef * x2
        x2 *= rf

    y = s2 - s1

    return y