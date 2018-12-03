from DadosExperimentais.LeituraDados import read_AST, read_LineDisplay, read_LineDisplayGroup
from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import time
import math as m
import matplotlib.pyplot as plt

# Leitura dos dados experimentais
# blocks = read_LineDisplay("Linha1_line_display", 0)
blocks = read_AST("AST_Samos")
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
Locate.SetVelocity(5.0)  # mm / us = km / s
# Para o AST v = 5km/s --> onda superficial
# Para quebra de grafite v = 3.2 km/s --> onda longitudinal


# Posição dos sensores
Locate.AddSensor(0, h / 3)  # 1
Locate.AddSensor(C / 3, h / 3)  # 2
Locate.AddSensor(2 * C / 3, h / 3)  # 3
Locate.AddSensor(C / 6, 2130.0)  # 4
Locate.AddSensor(C / 2, 2130.0)  # 5
Locate.AddSensor(5 * C / 6, 2130.0)  # 6
Locate.AddSensor(C / 2, -sp / 2)  # 7
Locate.AddSensor(C, -sp / 2)  # 8
Locate.AddSensor(0, h + sp / 2)  # 9
Locate.AddSensor(C / 2, h + sp / 2)  # 10

# Ponto de teste
k = 0
x_real = []
y_real = []
x_calc = []
y_calc = []
x_simple = []
y_simple = []
x_error = []
y_error = []
x_errorS = []
y_errorS = []

for block in blocks:
    # xp, yp = (block.X, block.Y)
    xp, yp = Locate.GetSensorCoords(block.pulser)
    x_real.append(xp)
    y_real.append(yp)
    k += 1

    read_data = zip(block.IDs, block.dt)
    data = []
    for (ID, dt) in read_data:
        data.append((ID, dt))

    print("Real: x: " + str(round(xp, 4)) + " / y: " + str(round(yp, 4)))

    x_s = Locate.simpleLocation(data)
    x_c = Locate.completeLocation(data)
    xs, ys = x_s[0], x_s[1]
    xc, yc = x_c[0], x_c[1]

    e_s = Locate.ExternalCalcDist(x_s[0], x_s[1], xp, yp)
    e_c = Locate.ExternalCalcDist(x_c[0], x_c[1], xp, yp)
    print("Erro do simplificado " + str(round(e_s, 3)) + " mm")
    print("Erro do completo " + str(round(e_c, 3)) + " mm")

    file.write('Ponto ' + str(k) + '\n')
    file.write("Real: x: " + str(round(xp, 4)) +
               " / y: " + str(round(yp, 4)) + '\n')
    file.write("Simples: x: " + str(round(xs, 4)) +
               " / y: " + str(round(ys, 4)) + '\n')
    file.write("Calculado: x: " + str(round(xc, 4)) +
               " / y: " + str(round(yc, 4)) + '\n')
    file.write('\n')

    if e_s < 300 or e_c < 300:
        x_simple.append(x_s[0])
        y_simple.append(x_s[1])

        x_errorS += [xp, x_s[0], np.nan]
        y_errorS += [yp, x_s[1], np.nan]

        print("Simples: x: " + str(round(x_s[0], 4)
                                   ) + " / y: " + str(round(x_s[1], 4)))

        x_calc.append(x_c[0])
        y_calc.append(x_c[1])

        x_error += [xp, x_c[0], np.nan]
        y_error += [yp, x_c[1], np.nan]

        print("Calculado: x: " +
              str(round(x_c[0], 4)) + " / y: " + str(round(x_c[1], 4)))

    print("\n")

file.close()
x_vessel = [0, C, C, 0, 0, 0, C, C, 0, 0, C, C, 0, 0]
y_vessel = [0, 0, h, h, 0, -sp, -sp, 0, 0, h, h, h + sp, h + sp, h]
plt.plot(x_vessel, y_vessel, 'k')
plt.plot(x_real, y_real, 'g.', markersize=12)
plt.plot(x_calc, y_calc, 'b.', markersize=12)
plt.plot(x_simple, y_simple, 'r.', markersize=12)
plt.plot(x_error, y_error, 'k--', linewidth=1)
plt.plot(x_errorS, y_errorS, 'k--', linewidth=1)
plt.legend(['Vaso', 'Real', 'Calc', 'Simples'])
plt.show()
