from DadosExperimentais.LeituraDados import read_LineDisplayGroup
from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import time
import math as m
import matplotlib.pyplot as plt

# Leitura dos dados experimentais
blocks = read_LineDisplayGroup("linha2", 150.0)
file = open('temp.txt', 'w')

# Parâmetros do vaso
C = 2492.0
d = C / m.pi
h = 2700.0
sp = 490.0

# Configurações do algoritmo
Locate = CylindricalLocation(d, h)
Locate.set_semiPerimeter(sp)
Locate.setCalcMode('section')
Locate.setSectionMode('reg')
Locate.SetVelocity(3.2)  # mm / us = km / s


# Posição dos sensores
Locate.AddSensor(0, h / 3)  # 1
Locate.AddSensor(C / 3, h / 3)  # 2
Locate.AddSensor(2 * C / 3, h / 3)  # 3
Locate.AddSensor(C / 6, 2130.0)  # 4
Locate.AddSensor(C / 2, 2130.0)  # 5
Locate.AddSensor(5 * C / 6, 2130.0)  # 6
Locate.AddSensor(C / 2, -sp * 0.4)  # 7
Locate.AddSensor(C, -sp * 0.4)  # 8
Locate.AddSensor(0, h + sp / 2)  # 9
Locate.AddSensor(C / 2, h + sp / 2)  # 10

# Ponto de teste
k = 0
x_real = []
y_real = []
x_calc = []
y_calc = []
x_disp = []
y_disp = []
x_error = []
y_error = []
x_errorD = []
y_errorD = []

for block in blocks:
    xp, yp = (block.X, block.Y)
    xd, yd = (block.x_disp, block.y_disp)
    x_real.append(xp)
    y_real.append(yp)
    k += 1

    read_data = zip(block.IDs, block.dt)
    data = []
    for (ID, dt) in read_data:
        data.append((ID, dt))

    print("Real: x: " + str(round(xp, 4)) + " / y: " + str(round(yp, 4)))
    print("Disp: x: " + str(round(xd, 4)) + " / y: " + str(round(yd, 4)))

    x_c = Locate.completeLocation(data)
    print("Calculado: x: " +
          str(round(x_c[0], 4)) + " / y: " + str(round(x_c[1], 4)))

    e_d = Locate.ExternalCalcDist(xd, yd, xp, yp)
    e_c = Locate.ExternalCalcDist(x_c[0], x_c[1], xp, yp)
    print("Erro do Disp " + str(round(e_d, 3)) + " mm")
    print("Erro do calculado " + str(round(e_c, 3)) + " mm")

    file.write('Ponto ' + str(k) + '\n')
    file.write("Real: x: " + str(round(xp, 4)) +
               " / y: " + str(round(yp, 4)) + '\n')
    file.write("Disp: x: " + str(round(xd, 4)) +
               " / y: " + str(round(yd, 4)) + '\n')
    file.write("Calculado: x: " +
               str(round(x_c[0], 4)) + " / y: " + str(round(x_c[1], 4)) + '\n')
    file.write('\n')

    if e_d < 300 and e_c < 300:
        x_disp.append(xd)
        y_disp.append(yd)

        x_errorD += [xp, xd, np.nan]
        y_errorD += [yp, yd, np.nan]

        x_calc.append(x_c[0])
        y_calc.append(x_c[1])

        x_error += [xp, x_c[0], np.nan]
        y_error += [yp, x_c[1], np.nan]

    print("\n")

file.close()
x_vessel = [0, C, C, 0, 0, 0, C, C, 0, 0, C, C, 0, 0]
y_vessel = [0, 0, h, h, 0, -sp, -sp, 0, 0, h, h, h + sp, h + sp, h]
plt.plot(x_vessel, y_vessel, 'k')
plt.plot(x_real, y_real, '.', markersize=12)
plt.plot(x_calc, y_calc, '.', markersize=12)
plt.plot(x_disp, y_disp, '.', markersize=12)
plt.plot(x_error, y_error, 'k--', linewidth=1)
plt.plot(x_errorD, y_errorD, 'k--', linewidth=1)
plt.legend(['Vaso', 'Real', 'Calc', 'Disp'])
plt.show()
