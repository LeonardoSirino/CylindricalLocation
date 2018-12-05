from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import matplotlib.pyplot as plt
import copy

file = open('Resultados\\AST_Samos.txt', 'r')

# Parâmetros do vaso
C = 2492.0
d = C / np.pi
h = 2700.0
sp = 490.0
Diag = np.sqrt((C / 2)**2 + h**2)

# Configurações do algoritmo
Locate = CylindricalLocation(d, h)
Locate.set_semiPerimeter(sp)
f = Locate.f
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
x_simple = []
y_simple = []

line = " "
while line != "":
    line = file.readline()
    if line[:4] == 'Real':
        ki = line.find(' / y: ')
        x = float(line[9:ki])
        y = float(line[ki + 6:])
        x_real.append(x)
        y_real.append(y)

    elif line[:7] == 'Simples':
        ki = line.find(' / y: ')
        x = float(line[12:ki])
        y = float(line[ki + 6:])
        x_simple.append(x)
        y_simple.append(y)

    elif line[:9] == 'Calculado':
        ki = line.find(' / y: ')
        x = float(line[14:ki])
        y = float(line[ki + 6:])
        x_calc.append(x)
        y_calc.append(y)

p_xreal = []
p_yreal = []
p_xsimple = []
p_ysimple = []
p_xcalc = []
p_ycalc = []
x_error = []
y_error = []
x_errorD = []
y_errorD = []
error_c = []
error_s = []

k = 0
for xp, yp, xs, ys, xc, yc in zip(x_real, y_real, x_simple, y_simple, x_calc, y_calc):
    k += 1
    e_s = Locate.ExternalCalcDist(xs, ys, xp, yp)
    e_c = Locate.ExternalCalcDist(xc, yc, xp, yp)
    error_c.append(e_c * 100 / Diag)
    error_s.append(e_s * 100 / Diag)
    print(str(k) + ' / ' + 'real: ' + str(xp) + ' / ' + str(yp) + ' / ' + 'simples: ' + str(xs) + ' / ' +
          str(ys) + ' / ' + str(e_s / Diag) + ' / ' + 'secc: ' + str(xc) + ' / ' + str(yc) + ' / ' + str(e_c / Diag))
    """
    print('Ponto ' + str(k))
    print("Erro do simples " + str(round(e_s, 3)) + " mm")
    print("Erro do calculado " + str(round(e_c, 3)) + " mm")
    """

    cond1 = e_s < 300 and e_c < 300
    cond1 = True
    cond2 = e_s < 0 or e_c < 0
    if cond1 or cond2:
        p_xreal.append(xp)
        p_yreal.append(yp)

        p_xsimple.append(xs)
        p_ysimple.append(ys)

        p_xcalc.append(xc)
        p_ycalc.append(yc)

        x_errorD += [xp, xs, np.nan]
        y_errorD += [yp, ys, np.nan]

        x_error += [xp, xc, np.nan]
        y_error += [yp, yc, np.nan]

    # print("\n")

file.close()


x_vessel = [0, C, C, 0, 0, 0, C, C, 0, 0, C, C, 0, 0]
y_vessel = [0, 0, h, h, 0, -sp, -sp, 0, 0, h, h, h + sp, h + sp, h]
plt.plot(x_vessel, y_vessel, 'k')
# plt.plot(x_sensor, y_sensor, 'y.', markersize=12)
plt.plot(x_real, y_real, 'g.', markersize=12)
plt.plot(p_xcalc, p_ycalc, 'b.', markersize=12)
plt.plot(p_xsimple, p_ysimple, 'r.', markersize=12)
plt.plot(x_error, y_error, 'b--', linewidth=1)
plt.plot(x_errorD, y_errorD, 'r--', linewidth=1)
plt.legend(['Vaso', 'Real', 'Seccionamento', 'Planificado'], loc=1)
plt.xlabel('Posição x [mm]')
plt.ylabel('Posição y [mm]')
plt.show()

plt.plot(y_real, error_s, 'r', y_real, error_c, 'b')
plt.legend(["Planificado", "Seccionamento"], loc=1)
plt.ylabel("Erro de posição [% diagonal do vaso]")
plt.xlabel("Altura [mm]")
plt.show()

