import geographiclib.geodesic as geo
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math as m
import numpy as np


class CylindricalLocation():
    def __init__(self, diameter, height, SemiPerimeter):
        self.diameter = diameter
        self.height = height
        self.SemiPerimeter = SemiPerimeter
        self.f = SemiPerimeter * 4 * 2 / \
            (m.pi * diameter) - 1  # Achatamento do tampo
        self.SensorList = []
        self.__SensorListHClone = []
        self.cap = geo.Geodesic(self.diameter / 2, self.f)

    def __AuxCoords(self, Coords):
        """
        Função para calcular as coordenadas auxliares (latitude e longitude) quando a coordenada estiver no tampo
        """
        if Coords.get("Ycord") >= self.height or Coords.get("Ycord") <= 0:
            if Coords.get("Xcord") > self.diameter * m.pi:
                AuxXcord = Coords.get("Xcord") - self.diameter * m.pi
            else:
                AuxXcord = Coords.get("Xcord")
            lon = AuxXcord / (self.diameter * m.pi) * 360 - 180
            if Coords.get("Ycord") >= self.height:
                Coords.update({"cap": "sup"})
                s12 = Coords.get("Ycord") - self.height
            else:
                s12 = abs(Coords.get("Ycord"))
                Coords.update({"cap": "inf"})

            EndCap = self.cap.Direct(lat1=0, lon1=lon, s12=s12, azi1=0)
            lat = EndCap.get("lat2")
            Coords.update({"lat": lat, "lon": lon, "OnCap": True})
        else:
            Coords.update({"OnCap": False})

    def __calcDist(self, P1, P2):
        """
        Função para calcular distância entre dois pontos em posições quaisquer do vaso
        """
        if P1.get("OnCap") and P2.get("OnCap"):
            if P1.get("cap") == P2.get("cap"):  # Dois pontos no mesmo tampo
                res = self.cap.Inverse(lat1=P1.get("lat"), lat2=P2.get(
                    "lat"), lon1=P1.get("lon"), lon2=P2.get("lon"))
                dist = res.get("s12")
            else:  # Pontos em tampos opostos
                # Definir os pontos corretos!
                dist = self.__DistCaptoCap(P1, P2)
        # Distância entre um ponto no casco e outro no tampo
        elif P1.get("OnCap") ^ P2.get("OnCap"):
            dist = self.__DistWalltoCap(P1, P2)
        else:  # Distãncia entre pontos no casco
            dist = m.sqrt((P1.get("Xcord") - P2.get("Xcord")) **
                          2 + (P1.get("Ycord") - P2.get("Ycord"))**2)

        return dist

    def __DistWalltoCap(self, P1, P2):
        if P1.get("OnCap"):
            Pcap = P1
            Pwall = P2
        else:
            Pcap = P2
            Pwall = P1

        if Pcap.get("cap") == "sup":
            Yaux = self.height
        else:
            Yaux = 0

        min = np.min([Pcap.get("Xcord"), Pwall.get("Xcord")])
        max = np.max([Pcap.get("Xcord"), Pwall.get("Xcord")])
        AuxPoint = {"Ycord": Yaux}
        distances = []

        for Xaux in np.linspace(start=min, stop=max, num=100):
            AuxPoint.update({"Xcord": Xaux})
            self.__AuxCoords(AuxPoint)
            AuxPoint.update({"OnCap": False})
            dist1 = self.__calcDist(Pwall, AuxPoint)
            AuxPoint.update({"OnCap": True})
            dist2 = self.__calcDist(AuxPoint, Pcap)
            distances.append(dist1 + dist2)

        dist = np.min(distances)

        return dist

    def __DistCaptoCap(self, Psup, Pinf):  # PlaceHolder
        print("Pontos em tampos opostos")
        dist = 0
        return dist

    def __DistVClone(self, Source, Sensor):  # PlaceHolder
        dist = 0
        return dist

    def AddSensor(self, Xcord, Ycord):
        # Conditions
        C1 = Xcord > 0 and Xcord < self.diameter * m.pi
        C2 = Ycord > - self.SemiPerimeter * \
            1.01 and Ycord < (self.height + self.SemiPerimeter) * 1.01
        if C1 and C2:
            SensorCoords = {'Xcord': Xcord, 'Ycord': Ycord}
            self.__AuxCoords(SensorCoords)
            self.SensorList.append(SensorCoords)
            if SensorCoords.get("Xcord") >= self.diameter * m.pi / 2:
                XHClone = SensorCoords.get("Xcord") - self.diameter * m.pi
            else:
                XHClone = SensorCoords.get("Xcord") + self.diameter * m.pi
            SensorHCloneCoords = {'Xcord': XHClone, 'Ycord': Ycord}
            self.__AuxCoords(SensorHCloneCoords)
            self.__SensorListHClone.append(SensorHCloneCoords)
            # Determining XVClone coordinate
            if Xcord > self.diameter / 2:
                Xcord = Xcord - self.diameter * m.pi / 2
            else:
                Xcord = Xcord + self.diameter * m.pi / 2
        else:
            print("As coordenadas deste ponto estão fora do vaso")

    def calcAllDist(self, SourceX, SourceY):
        i = 0
        Source = {"Xcord": SourceX, "Ycord": SourceY}
        self.__AuxCoords(Source)
        for sensor in self.SensorList:
            dist = self.__calcDist(Source, sensor)
            print(dist)
