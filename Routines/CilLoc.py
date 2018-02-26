import geographiclib.geodesic as geo
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math as m
import numpy as np


class VesselPoint():
    """Classe para definir as propriedades de um ponto no vaso
    """

    def __init__(self, Xcord, Ycord):
        self.Xcord = Xcord
        self.Ycord = Ycord
        self.Valid = True

    def SetXcord(self, Xcord):
        self.Xcord = Xcord

    def SetLon(self, Lon):
        self.Lon = Lon

    def SetLat(self, Lat):
        self.Lat = Lat

    def SetOnCap(self, OnCap):
        self.OnCap = OnCap

    def SetCap(self, Cap):
        self.Cap = Cap

    def SetValid(self, Valid):
        """Validade de um ponto para se calcular a distância.
        Os clones horizontais de pontos no tampo não são válidos e portanto não se deve a calcular a ditância até eles.
        Esses são mantidos para manter a equivalência entre os vetores de sensores e clones

        Arguments:
            Valid {[boolean]} -- [validade do ponto]
        """

        self.Valid = Valid


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
        if Coords.Ycord >= self.height or Coords.Ycord <= 0:
            if Coords.Xcord > self.diameter * m.pi:
                AuxXcord = Coords.Xcord - self.diameter * m.pi
            else:
                AuxXcord = Coords.Xcord
            lon = AuxXcord / (self.diameter * m.pi) * 360 - 180
            if Coords.Ycord >= self.height:
                Coords.SetCap("sup")
                s12 = Coords.Ycord - self.height
            else:
                s12 = abs(Coords.Ycord)
                Coords.SetCap("inf")

            EndCap = self.cap.Direct(lat1=0, lon1=lon, s12=s12, azi1=0)
            lat = EndCap.get("lat2")
            Coords.SetLat(lat)
            Coords.SetLon(lon)
            Coords.SetOnCap(True)
        else:
            Coords.SetOnCap(False)

    def __calcDist(self, P1, P2):
        """
        Função para calcular distância entre dois pontos em posições quaisquer do vaso
        """
        if P1.OnCap and P2.OnCap:
            if P1.Cap == P2.Cap:  # Dois pontos no mesmo tampo
                res = self.cap.Inverse(
                    lat1=P1.Lat, lat2=P2.Lat, lon1=P1.Lon, lon2=P2.Lon)
                dist = res.get("s12")
            else:  # Pontos em tampos opostos
                # Definir os pontos corretos!
                dist = self.__DistCaptoCap(P1, P2)
        # Distância entre um ponto no casco e outro no tampo
        elif P1.OnCap ^ P2.OnCap:
            dist = self.__DistWalltoCap(P1, P2)
        else:  # Distãncia entre pontos no casco
            dist = m.sqrt((P1.Xcord - P2.Xcord) **
                          2 + (P1.Ycord - P2.Ycord)**2)

        return dist

    def __DistWalltoCap(self, P1, P2):
        if P1.OnCap:
            Pcap = P1
            Pwall = P2
        else:
            Pcap = P2
            Pwall = P1

        if Pcap.Cap == "sup":
            Yaux = self.height
        else:
            Yaux = 0

        min = np.min([Pcap.Xcord, Pwall.Xcord])
        max = np.max([Pcap.Xcord, Pwall.Xcord])
        AuxPoint = VesselPoint(0, Yaux)
        distances = []

        for Xaux in np.linspace(start=min, stop=max, num=100):
            AuxPoint.SetXcord(Xaux)
            self.__AuxCoords(AuxPoint)
            AuxPoint.SetOnCap(False)
            dist1 = self.__calcDist(Pwall, AuxPoint)
            AuxPoint.SetOnCap(True)
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
            SensorCoords = VesselPoint(Xcord, Ycord)
            self.__AuxCoords(SensorCoords)
            self.SensorList.append(SensorCoords)
            SensorHCloneCoords = VesselPoint(0, Ycord)
            if SensorCoords.Ycord >= 0 and SensorCoords.Ycord <= self.height:  # Ponto no casco do vaso
                # Sensor posicionado na porção direita do casco
                if SensorCoords.Xcord >= self.diameter * m.pi / 2:
                    XHClone = SensorCoords.Xcord - self.diameter * m.pi
                else:
                    # Sensor posicionado na porção esquerda do casco
                    XHClone = SensorCoords.Xcord + self.diameter * m.pi
                    SensorHCloneCoords.SetXcord(XHClone)
            else:  # Ponto em um dos tampos do vaso
                SensorHCloneCoords.SetXcord(m.nan)
                SensorHCloneCoords.SetValid(False)

            self.__AuxCoords(SensorHCloneCoords)
            self.__SensorListHClone.append(SensorHCloneCoords)
        else:
            print("As coordenadas deste ponto estão fora do vaso")

    def calcAllDist(self, SourceX, SourceY):
        Source = VesselPoint(SourceX, SourceY)
        self.__AuxCoords(Source)
        for sensor in self.SensorList:
            dist = self.__calcDist(Source, sensor)
            print(dist)
