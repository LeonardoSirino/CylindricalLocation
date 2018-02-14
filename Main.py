import geographiclib.geodesic as geo
import matplotlib.pyplot as plt
import math as m


class CylindricalLocation():
    def __init__(self, diameter, height, SemiPerimeter):
        self.diameter = diameter
        self.height = height
        self.SemiPerimeter = SemiPerimeter
        self.SensorList = []
        self.SensorListHClone = []
        self.SensorListVClone = []

    def AddSensor(self, Xcord, Ycord):
        # Conditions
        C1 = Xcord > 0 and Xcord < self.diameter
        C2 = Ycord > - \
            self.diameter and Ycord < (self.height + self.SemiPerimeter)
        if C1 and C2:
            self.SensorList.append({'Xcord': Xcord, 'Ycord': Ycord})
            self.SensorListHClone.append(
                {'Xcord': Xcord + self.diameter, 'Ycord': Ycord})
            # Determining XVClone coordinate
            if Xcord > self.diameter / 2:
                Xcord = Xcord - self.diameter / 2
            else:
                Xcord = Xcord + self.diameter / 2
            if Ycord > self.height / 2:
                self.SensorListVClone.append(
                    {'Xcord': Xcord, 'Ycord': Ycord + 2 * self.SemiPerimeter + (self.height - Ycord)})
            else:
                self.SensorListVClone.append(
                    {'Xcord': Xcord, 'Ycord': Ycord - 2 * self.SemiPerimeter - (self.height - Ycord)})
        else:
            print("As coordenadas deste ponto estão fora do vaso")

    def calcAllDist(self, SourceX, SourceY):
        i=0
        for Original in self.SensorList:
            print(Original)
            HClone=self.SensorListHClone[i]
            VClone=self.SensorListVClone[i]
            distOriginal = m.sqrt(
                (Original.get('Xcord') - SourceX) ** 2 + (Original.get('Ycord') - SourceY) ** 2)
            print(distOriginal)
            distHClone = m.sqrt(
                (HClone.get('Xcord') - SourceX) ** 2 + (HClone.get('Ycord') - SourceY) ** 2)
            print(distHClone)
            distVClone = m.sqrt(
                (VClone.get('Xcord') - SourceX) ** 2 + (VClone.get('Ycord') - SourceY) ** 2)
            print(distVClone)
            i+=1

Locate=CylindricalLocation(1000.0,1000.0,100.0)
Locate.AddSensor(Xcord=200.0,Ycord=200.0)
Locate.calcAllDist(100.0,100.0)