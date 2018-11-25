from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import time
import math as m
import matplotlib.pyplot as plt

# Parâmetros do vaso
C = 2492.0
diameter = C / m.pi
height = 2700.0
semiperimeter = 510.0

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
Locate.AddSensor(C / 6, 2 * h / 3)  # 4
Locate.AddSensor(C / 2, 2 * h / 3)  # 5
Locate.AddSensor(5 * C / 6, 2 * h / 3)  # 6
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
divs = 5
x_array = np.linspace(diameter * m.pi * 0.05, diameter * m.pi * 0.95, num=divs)

# Inicialização dos vetores
x_RP = np.zeros(divs)
y_RP = np.zeros(divs)
x_IK = []
y_IK = []
x_SL = []
y_SL = []
x_CL = []
y_CL = []
error_IK = []
error_S = []
error_Sec = []

j = 0
yp = 0
for xp in x_array:
    t = Locate.returnDeltaT(xp, yp, [-1], 'geodesic')
    data = []
    i = 0
    for AT in t:
        data.append((i, AT))
        i += 1

    x_RP[j] = xp
    y_RP[j] = yp

    x = Locate._CylindricalLocation__InitialKick(data)
    x_IK += [xp, x[0], np.nan]
    y_IK += [yp, x[1], np.nan]
    e = Locate.ExternalCalcDist(xp, yp, x[0], x[1])
    error_IK.append(e)

    x = Locate.simpleLocation(data)
    x_SL += [xp, x[0], np.nan]
    y_SL += [yp, x[1], np.nan]
    e = Locate.ExternalCalcDist(xp, yp, x[0], x[1])
    error_S.append(e)

    x = Locate.completeLocation(data)
    x_CL += [xp, x[0], np.nan]
    y_CL += [yp, x[1], np.nan]
    e = Locate.ExternalCalcDist(xp, yp, x[0], x[1])
    error_Sec.append(e)

    j += 1
    print(j)

x_vessel = [0, C, C, 0, 0, 0, C, C, 0, 0, C, C, 0, 0]
y_vessel = [0, 0, h, h, 0, -sp, -sp, 0, 0, h, h, h + sp, h + sp, h]
plt.plot(x_vessel, y_vessel, 'k')
plt.plot(x_RP, y_RP, 'go', x_SL, y_SL, x_CL, y_CL, xS, yS, 'yo')
plt.legend(["Vaso", "Posição real", "Simples", "Completa", "Sensores"])
plt.xlabel("Coordenada X")
plt.ylabel("Coordenada Y")
plt.show()

plt.plot(x_RP, error_S, x_RP, error_Sec)
plt.legend(["Simples", "Completa"])
plt.ylabel("Erro de posição [mm]")
plt.xlabel("Posição na interface [mm]")
plt.show()