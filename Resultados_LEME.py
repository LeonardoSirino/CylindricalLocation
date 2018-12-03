from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import matplotlib.pyplot as plt

file = open('Resultados\\linha 3.txt', 'r')

# Parâmetros do vaso
C = 2492.0
d = C / np.pi
h = 2700.0
sp = 490.0

# Configurações do algoritmo
Locate = CylindricalLocation(d, h)
Locate.set_semiPerimeter(sp)
Locate.setCalcMode('section')
Locate.setSectionMode('reg')
Locate.SetVelocity(3.2)  # mm / us = km / s

# Posição dos sensores
x_sensor = [0, C / 3, 2 * C / 3, C / 6, C / 2, 5 * C / 6, C / 2, C, 0, C / 2]
y_sensor = [h / 3, h / 3, h / 3, 2130.0, 2130.0,
            2130.0, -sp * 0.5, -sp * 0.5, h + sp / 2, h + sp / 2]

# Ponto de teste
k = 0
x_real = []
y_real = []
x_calc = []
y_calc = []
x_disp = []
y_disp = []

line = " "
while line != "":
    line = file.readline()
    if line[:4] == 'Real':
        ki = line.find(' / y: ')
        x = float(line[9:ki])
        y = float(line[ki + 6:])
        x_real.append(x)
        y_real.append(y)

    elif line[:4] == 'Disp':
        ki = line.find(' / y: ')
        x = float(line[9:ki])
        y = float(line[ki + 6:])
        x_disp.append(x)
        y_disp.append(y)

    elif line[:9] == 'Calculado':
        ki = line.find(' / y: ')
        x = float(line[14:ki])
        y = float(line[ki + 6:])
        x_calc.append(x)
        y_calc.append(y)

p_xreal = []
p_yreal = []
p_xdisp = []
p_ydisp = []
p_xcalc = []
p_ycalc = []
x_error = []
y_error = []
x_errorD = []
y_errorD = []
error_c = []
error_d = []
k = 0
for xp, yp, xd, yd, xc, yc in zip(x_real, y_real, x_disp, y_disp, x_calc, y_calc):
    k += 1
    e_d = Locate.ExternalCalcDist(xd, yd, xp, yp)
    e_c = Locate.ExternalCalcDist(xc, yc, xp, yp)
    print('Ponto ' + str(k))
    print("Erro do Disp " + str(round(e_d, 3)) + " mm")
    print("Erro do calculado " + str(round(e_c, 3)) + " mm")

    cond1 = e_d < 300 and e_c < 300
    cond2 = e_d < 0 or e_c < 0
    if cond1 or cond2:
        p_xreal.append(xp)
        p_yreal.append(yp)

        p_xdisp.append(xd)
        p_ydisp.append(yd)

        p_xcalc.append(xc)
        p_ycalc.append(yc)

        x_errorD += [xp, xd, np.nan]
        y_errorD += [yp, yd, np.nan]

        x_error += [xp, xc, np.nan]
        y_error += [yp, yc, np.nan]

        error_c.append(e_c)
        error_d.append(e_d)

    print("\n")

file.close()

x_vessel = [0, C, C, 0, 0, 0, C, C, 0, 0, C, C, 0, 0]
y_vessel = [0, 0, h, h, 0, -sp, -sp, 0, 0, h, h, h + sp, h + sp, h]
plt.plot(x_vessel, y_vessel, 'k')
plt.plot(x_sensor, y_sensor, '.', markersize=12)
plt.plot(x_real, y_real, '.', markersize=12)
plt.plot(p_xcalc, p_ycalc, '.', markersize=12)
plt.plot(p_xdisp, p_ydisp, '.', markersize=12)
plt.plot(x_error, y_error, 'k--', linewidth=1)
plt.plot(x_errorD, y_errorD, 'k--', linewidth=1)
plt.legend(['Vaso', 'Sensores', 'Real', 'Calc', 'Disp'])
plt.show()

plt.plot(p_xreal, error_d, '.', p_xreal, error_c, '.')
plt.legend(['Disp', 'Calculado'])
plt.ylim((0, 300))
plt.show()
