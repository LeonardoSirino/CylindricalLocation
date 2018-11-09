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

@jit(nopython=True, parallel=True)
def distCap(redF, x1, y1, x2, y2, tol, d):
	r1q = x1**2 + y1**2
	if r1q - d**2 >= 0:
		u1 = np.sqrt(r1q - d**2)
	else:
		if abs(r1q - d) < tol:
			u1 = 0
		else:
			u1 = np.nan

	r2q = x2**2 + y2**2
	if r2q - d**2 >= 0:
		u2 = np.sqrt(r2q - d**2)
	else:
		if abs(r2q - d) < tol:
			u2 = 0
		else:
			u2 = np.nan

	v1c = np.array([-x1, -y1])
	v2c = np.array([-x2, -y2])
	v12 = np.array([x2 - x1, y2 - y1])

	Nv1c = np.linalg.norm(v1c)
	Nv2c = np.linalg.norm(v2c)
	Nv12 = np.linalg.norm(v12)

	if Nv1c > 0 and Nv12 > 0:
		cosTheta1 = np.dot(v1c, v12) / (Nv1c * Nv12)
	else:
		cosTheta1 = 0

	if Nv2c > 0 and Nv12 > 0:
		cosTheta2 = np.dot(v2c, v12) / (Nv2c * Nv12)
	else:
		cosTheta2 = 0

	if cosTheta1 <= 0:
		# Ponto está do outro lado do centro da elipse
		u1 = -u1

	if cosTheta2 <= 0:
		# Ponto está do outro lado do centro da elipse
		u2 = -u2

	return u1, u2