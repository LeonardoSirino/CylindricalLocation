from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import time
import math as m
import matplotlib.pyplot as plt

# Parâmetros do vaso
diameter = 400.0
height = 1000.0
f = 0.5

# Configurações do algoritmo
Locate = CylindricalLocation(diameter, height)
Locate.setCalcMode('section')
Locate.set_f(f)
semiperimeter = Locate.SemiPerimeter
Locate.SetVelocity(5)

# Posição dos sensores
Locate.StructuredSensorDistribution(lines=3, sensorsInLine=3, x0=0, y0=height +
                                    semiperimeter / 4, dx=(diameter * m.pi) / 4, dy=semiperimeter / 4, aligned=False)

# Vetores de posição dos sensores
xS = []
yS = []
for sensor in Locate.SensorList:
    xS.append(sensor.Xcord)
    yS.append(sensor.Ycord)

# Pontos de teste
xdivs = 5
ydivs = 5
x_array = np.linspace(diameter * m.pi * 0.01, diameter * m.pi * 0.99, num=xdivs)
y_array = np.linspace(height * 1.1, height + semiperimeter, num=ydivs)

# Inicialização dos vetores
x_RP = np.zeros(xdivs * ydivs)
y_RP = np.zeros(xdivs * ydivs)
x_IK = []
y_IK = []
x_SL = []
y_SL = []
x_CL = []
y_CL = []

j = 0
for xp in x_array:
    for yp in y_array:
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

        x = Locate.simpleLocation(data)
        x_SL += [xp, x[0], np.nan]
        y_SL += [yp, x[1], np.nan]

        x = Locate.completeLocation(data)
        x_CL += [xp, x[0], np.nan]
        y_CL += [yp, x[1], np.nan]

        j += 1
        print(j)

plt.plot(x_RP, y_RP, 'ko', x_IK, y_IK, x_SL, y_SL, x_CL, y_CL, xS, yS, 'yo')
plt.legend(["Posição real", "Chute inicial", "Simples", "Completa", "Sensores"])
plt.xlabel("Coordenada X")
plt.ylabel("Coordenada Y")
plt.ylim((900, 1500))
plt.xlim((-100, 1400))
plt.show()