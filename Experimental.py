from DadosExperimentais.LeituraDados import read_AST
from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import time
import math as m

# Leitura dos dados experimentais
blocks = read_AST("AST")

# Parâmetros do vaso
C = 2492.0
d = C / m.pi
h = 2700.0
f = 0.5

# Configurações do algoritmo
Locate = CylindricalLocation(d, h)
Locate.setCalcMode('section')
Locate.set_f(f)
sp = Locate.SemiPerimeter
Locate.SetVelocity(4.880) # mm / uS

# Posição dos sensores
Locate.AddSensor(0, 0) # Não usado
Locate.AddSensor(0, h / 3) # 1
Locate.AddSensor(C / 3, h / 3) # 2
Locate.AddSensor(2 * C / 3, h / 3) # 3
Locate.AddSensor(C / 6, 2 * h / 3) # 4
Locate.AddSensor(C / 2, 2 * h / 3) # 5
Locate.AddSensor(5 * C / 6, 2 * h / 3) # 6
Locate.AddSensor(C / 2, -sp / 2) # 7
Locate.AddSensor(C, -sp / 2) # 8
Locate.AddSensor(0, h + sp / 2) # 9
Locate.AddSensor(C / 2, h + sp / 2) # 10

# Ponto de teste
block = blocks[4]
xp = C / 6
yp = 2 * h / 3

t = block.dt
data = []
i = 0
for AT in t:
    i += 1
    data.append((i, AT))

print(data)

print("Posição verdadeira:")
print("x: " + str(round(xp, 4)) +
      " / y: " + str(round(yp, 4)))
print("\n")

t0 = time.time()

x = Locate.simpleLocation(data)

t1 = time.time()
print("Localização simplificada:")
print(x)
print("Tempo decorrido: " + str(t1 - t0))
print("\n")

t0 = time.time()

x = Locate.completeLocation(data)

t1 = time.time()
print("Localização completa:")
print(x)
print("Tempo decorrido: " + str(t1 - t0))
