from Routines.CilLoc import VesselPoint, CylindricalLocation
import numpy as np
import matplotlib.pyplot as plt
import copy

file = open('Resultados\\linha 3.txt', 'r')


def lenghtByHeight(z0):
    if z0 < 0:
        global d
        global f

        s = 0
        z = 0
        x = -d / 2
        a = d / 2
        dx = a * 0.00005
        while z < abs(z0):
            xi = x + dx
            zi = np.sqrt(1 - xi**2 / a**2) * a * f
            dz = zi - z
            ds = np.sqrt(dx**2 + dz**2)
            s += ds
            z = zi
            x = xi

        return -s

    else:
        return z0


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
x_array = np.linspace(0, 2250, num=10)
data_disp = [[]] * 10
data_calc = [[]] * 10

k = 0
for xp, yp, xd, yd, xc, yc in zip(x_real, y_real, x_disp, y_disp, x_calc, y_calc):
    yd = lenghtByHeight(yd)
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

        i = np.where(x_array == xp)[0][0]
        data = copy.copy(data_disp[i])
        data.append(e_d)
        data_disp[i] = data

        data = copy.copy(data_calc[i])
        data.append(e_c)
        data_calc[i] = data

    print("\n")

file.close()

x_vessel = [0, C, C, 0, 0, 0, C, C, 0, 0, C, C, 0, 0]
y_vessel = [0, 0, h, h, 0, -sp, -sp, 0, 0, h, h, h + sp, h + sp, h]
plt.plot(x_vessel, y_vessel, 'k')
plt.plot(x_sensor, y_sensor, 'y.', markersize=12)
plt.plot(x_real, y_real, 'g.', markersize=12)
plt.plot(p_xcalc, p_ycalc, 'b.', markersize=12)
plt.plot(p_xdisp, p_ydisp, 'r.', markersize=12)
plt.plot(x_error, y_error, 'k--', linewidth=1)
plt.plot(x_errorD, y_errorD, 'k--', linewidth=1)
plt.legend(['Vaso', 'Sensores', 'Real', 'Seccionamento', 'AEWin'], loc=1)
plt.xlabel('Posição x [mm]')
plt.ylabel('Posição y [mm]')
plt.show()

x_vessel = [0, C, C, 0, 0, 0, C, C, 0, 0, C, C, 0, 0]
y_vessel = [0, 0, h, h, 0, -sp, -sp, 0, 0, h, h, h + sp, h + sp, h]
plt.plot(x_vessel, y_vessel, 'k')
plt.plot(x_sensor, y_sensor, 'y.', markersize=12)
plt.plot(x_real, y_real, 'g.', markersize=12)
plt.plot(p_xcalc, p_ycalc, 'b.', markersize=12)
plt.plot(p_xdisp, p_ydisp, 'r.', markersize=12)
plt.plot(x_error, y_error, 'k--', linewidth=1)
plt.plot(x_errorD, y_errorD, 'k--', linewidth=1)
plt.legend(['Vaso', 'Sensores', 'Real', 'Seccionamento', 'AEWin'], loc=1)
plt.ylim((-300, 300))
plt.xlabel('Posição x [mm]')
plt.ylabel('Posição y [mm]')
plt.show()

mean_c = []
mean_d = []
dp_c = []
dp_d = []
for d_c, d_d in zip(data_calc, data_disp):
    mean_c.append(np.mean(d_c) * 100 / Diag)
    dp_c.append(np.std(d_c) * 100 / Diag)
    mean_d.append(np.mean(d_d) * 100 / Diag)
    dp_d.append(np.std(d_d) * 100 / Diag)

plt.errorbar(x_array - [10] * len(x_array), mean_d,
             yerr=dp_d, uplims=True, lolims=True, fmt='ro')
plt.errorbar(x_array + [10] * len(x_array), mean_c,
             yerr=dp_c, uplims=True, lolims=True, fmt='bo')
plt.legend(['AEWin', 'Seccionamento'], loc=1)
plt.ylabel('Erro [% diagonal]')
plt.xlabel('Posição X [mm]')
plt.show()
