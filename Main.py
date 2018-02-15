import geographiclib.geodesic as geo
import matplotlib.pyplot as plt
import math as m


class CylindricalLocation():
    def __init__(self, diameter, height, SemiPerimeter):
        self.diameter = diameter
        self.height = height
        self.SemiPerimeter = SemiPerimeter
        self.f = SemiPerimeter * 4 * 2 / \
            (m.pi * diameter) - 1  # Achatamento do tampo
        self.SensorList = []
        self.SensorListHClone = []
        self.SensorListVClone = []

    def AddSensor(self, Xcord, Ycord):
        # Conditions
        C1 = Xcord > 0 and Xcord < self.diameter * m.pi
        C2 = Ycord > - \
            self.SemiPerimeter and Ycord < (self.height + self.SemiPerimeter)
        if C1 and C2:
            self.SensorList.append({'Xcord': Xcord, 'Ycord': Ycord})
            self.SensorListHClone.append(
                {'Xcord': Xcord + self.diameter * m.pi, 'Ycord': Ycord})
            # Determining XVClone coordinate
            if Xcord > self.diameter / 2:
                Xcord = Xcord - self.diameter * m.pi / 2
            else:
                Xcord = Xcord + self.diameter * m.pi / 2
            # Determining YVClone coordinate
            if Ycord > self.height:  # No tampo superior
                Ycord = Ycord + (self.height + self.SemiPerimeter - Ycord)
            elif Ycord > self.height / 2 and Ycord < self.height:  # Na parte superior do corpo
                Ycord = Ycord + 2 * self.SemiPerimeter + (self.height - Ycord)
            elif Ycord < self.height / 2 and Ycord > 0:  # Na parte inferior do corpo
                Ycord = Ycord - 2 * self.SemiPerimeter - (self.height - Ycord)
            else:  # No tampo inferior
                Ycord = Ycord - (self.SemiPerimeter - abs(Ycord))
            # Setting VClone Coordinates
            self.SensorListVClone.append({'Xcord': Xcord, 'Ycord': Ycord})

        else:
            print("As coordenadas deste ponto estão fora do vaso")

    def calcAllDist(self, SourceX, SourceY):
        i = 0
        for Original in self.SensorList:
            HClone = self.SensorListHClone[i]
            VClone = self.SensorListVClone[i]
            distOriginal = m.sqrt(
                (Original.get('Xcord') - SourceX) ** 2 + (Original.get('Ycord') - SourceY) ** 2)
            print("Distância original: " + str(distOriginal))
            distHClone = m.sqrt(
                (HClone.get('Xcord') - SourceX) ** 2 + (HClone.get('Ycord') - SourceY) ** 2)
            print("Distância clone horizontal: " + str(distHClone))
            distVClone = m.sqrt(
                (VClone.get('Xcord') - SourceX) ** 2 + (VClone.get('Ycord') - SourceY) ** 2)
            print("Distância clone vertical: " + str(distVClone))
            i += 1


Locate = CylindricalLocation(1000.0, 1000.0, 750 * m.pi / 4)
Locate.AddSensor(Xcord=200.0, Ycord=200.0)
Locate.calcAllDist(100.0, 100.0)
