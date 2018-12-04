from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import time
import math as m
import matplotlib.pyplot as plt

file = open('temp.txt', 'w')

# Parâmetros do vaso
C = 2492.0
diameter = C / m.pi
height = 2700.0
semiperimeter = 490.0

h = height
sp = semiperimeter

# Configurações do algoritmo
Locate = CylindricalLocation(diameter, height)
Locate.setCalcMode('section')
Locate.set_semiPerimeter(semiperimeter)
Locate.SetVelocity(5)

# Posição dos sensores
Locate.AddSensor(0, h / 3)  # 1
Locate.AddSensor(C / 3, h / 3)  # 2
Locate.AddSensor(2 * C / 3, h / 3)  # 3
Locate.AddSensor(C / 6, 2130)  # 4
Locate.AddSensor(C / 2, 2130)  # 5
Locate.AddSensor(5 * C / 6, 2130)  # 6
Locate.AddSensor(C / 2, -sp / 2)  # 7
Locate.AddSensor(C, -sp / 2)  # 8
Locate.AddSensor(0, h + sp / 2)  # 9
Locate.AddSensor(C / 2, h + sp / 2)  # 10

# Vetores de posição dos sensores
xS = []
yS = []
for sensor in Locate.SensorList:
    xS.append(sensor.Xcord)
    yS.append(sensor.Ycord)

# Pontos de teste
divs = 20
var_array = np.linspace(-0.95 * semiperimeter, height / 2, num=divs)

# Inicialização dos vetores
x_RP = np.zeros(divs)
y_RP = np.zeros(divs)
x_SL = []
y_SL = []
x_CL = []
y_CL = []
x_eSL = []
y_eSL = []
x_eCL = []
y_eCL = []
error_S = np.zeros(divs)
error_Sec = np.zeros(divs)

j = 0
xp = C / 2
for yp in var_array:
    t = Locate.returnDeltaT(xp, yp, [-1], 'geodesic')
    data = []
    i = 0
    for AT in t:
        data.append((i, AT))
        i += 1

    x_RP[j] = xp
    y_RP[j] = yp

    file.write('Ponto ' + str(j) + '\n')
    file.write("Real: x: " + str(round(xp, 4)) +
                   " / y: " + str(round(yp, 4)) + '\n')

    x = Locate.simpleLocation(data)
    e_s = Locate.ExternalCalcDist(x[0], x[1], xp, yp)
    x_eSL += [xp, x[0], np.nan]
    y_eSL += [yp, x[1], np.nan]
    x_SL.append(x[0])
    y_SL.append(x[1])
    file.write("Simples: x: " +
                str(round(x[0], 4)) + " / y: " + str(round(x[1], 4)) + '\n')
    error_S[j] = e_s

    x = Locate.completeLocation(data)
    e_c = Locate.ExternalCalcDist(x[0], x[1], xp, yp)
    x_eCL += [xp, x[0], np.nan]
    y_eCL += [yp, x[1], np.nan]
    x_CL.append(x[0])
    y_CL.append(x[1])
    file.write("Calculado: x: " +
                str(round(x[0], 4)) + " / y: " + str(round(x[1], 4)) + '\n')
    error_Sec[j] = e_c

    file.write('\n')

    j += 1
    print(j)

file.close()

x_vessel = [0, C, C, 0, 0, 0, C, C, 0, 0, C, C, 0, 0]
y_vessel = [0, 0, h, h, 0, -sp, -sp, 0, 0, h, h, h + sp, h + sp, h]
x_sensor = [0, C / 3, 2 * C / 3, C / 6, C / 2, 5 * C / 6, C / 2, C, 0, C / 2]
y_sensor = [h / 3, h / 3, h / 3, 2130.0, 2130.0,
            2130.0, -sp * 0.5, -sp * 0.5, h + sp / 2, h + sp / 2]
plt.plot(x_vessel, y_vessel, 'k')
plt.plot(x_sensor, y_sensor, 'y.', markersize=12)
plt.plot(x_RP, y_RP, 'g.', markersize=12)
plt.plot(x_SL, y_SL, 'r.', markersize=12)
plt.plot(x_CL, y_CL, 'b.', markersize=12)
plt.plot(x_eSL, y_eSL, 'r--', linewidth=1)
plt.plot(x_eCL, y_eCL, 'b--', linewidth=1)
plt.legend(['Vaso', 'Sensores', 'Real', 'Planificado', 'Seccionamento'], loc=1)
plt.xlabel('Posição x [mm]')
plt.ylabel('Posição y [mm]')
plt.show()

plt.plot(y_RP, error_S, y_RP, error_Sec)
plt.legend(["Simples", "Completa"])
plt.ylabel("Erro de posição [mm]")
plt.xlabel("Posição na interface [mm]")
plt.show()