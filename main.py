
from LatLon import lat_lon as ll
from obspy import read_events
from json import load
from utils.synEq import generateTT
from utils.catalog2hypo import catalog2hypoellipse, station2hypoellipse, velocity2hypoellipse, generateHypoellipseDefaultFile
from obspy.geodetics.base import gps2dist_azimuth as gps
from numpy import array, max
from pathlib import Path
from string import ascii_letters as al
from random import choices
import os

# define main class
class main():
    def __init__(self):
        with open("par.json") as f:
            self.parDict = load(f)
        with open("hypoDefaults.json") as f:
            self.hypoellipseDefaultsDict = load(f)            
        self.velocityModel = self.readVelocityFile(self.parDict["velocityFile"])
        self.stations = self.readVelocityFile(self.parDict["stationFile"])
        self.zGridMax = int(max(self.velocityModel["Z"]) + 1.0)
        self.node_intervals = self.parDict["velocityGridIntervals"]
        self.decimationFactor = self.parDict["decimationFactor"]

    def readVelocityFile(self, stationFile):
        """
        Load velocity model from "STATION0.HYP" file.
        - Inputs:
        stationFile: full name of NORDIC station file,
        - Output:
        velocityModel: a dictionary containing velocity model.
        """
        msg = "+++ Parsing velocity model ..."
        emptyLines = 0
        velocityModel = {"Vp": [], "Z": [], "VpVs": 1.73}
        with open(stationFile) as f:
            for l in f:
                if not l.strip():
                    emptyLines += 1
                if emptyLines == 2 and l.strip():
                    Vp, Z = [float(x) for x in l.split()[:2]]
                    velocityModel["Vp"].append(Vp)
                    velocityModel["Z"].append(Z)
                if emptyLines == 3 and l.strip():
                    VpVs = float(l[16:20])
                    velocityModel["VpVs"] = VpVs
                    break
        return velocityModel

    def readStationFile(self, stationFile):
        """
        Read station information from "STATION0.HYP" file.
        - Inputs:
        stationFile: full name of NORDIC station file,
        - Output:
        stations: a dictionary containing stations information.
        """    
        msg = "+++ Parsing station information ..."
        emptyLines = 0
        stations = {}
        with open(stationFile) as f:
            for l in f:
                if not l.strip():
                    emptyLines += 1
                if emptyLines == 1 and l.strip():
                    code = l[:6].strip()
                    lat = ll.Latitude(degree=int(
                        l[6:8]), minute=float(l[8:13])).decimal_degree
                    lon = ll.Longitude(degree=int(
                        l[15:17]), minute=float(l[17:22])).decimal_degree
                    elv = float(l[23:27])
                    stations[code] = {"Lat":lat, "Lon":lon, "Elv":elv}
        return stations

    def readEventFile(self):
        catalog = read_events(self.parDict["eventsFile"])
        return catalog

    
    def computeNewDistances(self, eventLat, eventLon, picks, stations):
        stationsInPicks = [pick.waveform_id.station_code for pick in picks]
        stationLats = [stations[station]["Lat"] for station in stationsInPicks]
        stationLons = [stations[station]["Lon"] for station in stationsInPicks]
        distances = [gps(eventLat, eventLon, stationLat, stationLon)[0]*1e-3 for (stationLat, stationLon) in zip(stationLats, stationLons)]
        return [array([distance, 1.0, 0.0]) for distance in distances]

    def writeCatalogFile(self, catalog, id):
        catPath = Path(os.path.join("tmp"))
        catPath.mkdir(parents=True, exist_ok=True)
        outFile = os.path.join(catPath, "{id:6s}.pha".format(id=id))
        catalog2hypoellipse(catalog, outFile)
    
    def writeStationFile(self, stations, id):
        staPath = Path(os.path.join("tmp"))
        staPath.mkdir(parents=True, exist_ok=True)
        outFile = os.path.join(staPath, "{id:6s}.sta".format(id=id))
        station2hypoellipse(stations, outFile)        

    def writeVelocityFile(self, velocity, id):
        velPath = Path(os.path.join("tmp"))
        velPath.mkdir(parents=True, exist_ok=True)
        outFile = os.path.join(velPath, "{id:6s}.prm".format(id=id))
        velocity2hypoellipse(velocity, outFile)

    def writeHypoellipseConfigFile(self):
        defPath = Path(os.path.join("tmp"))
        defPath.mkdir(parents=True, exist_ok=True)
        outFile = os.path.join(defPath, "default.cfg")
        generateHypoellipseDefaultFile(self.hypoellipseDefaultsDict, outFile)

    def computeNewArrivals(self, catalog, updatedStations):
        for event in catalog:
            preferred_origin = event.preferred_origin()
            originTime = preferred_origin.time
            eventLatitude = preferred_origin.latitude
            eventLongitude = preferred_origin.longitude
            depth = preferred_origin.depth
            picksP = [pick for pick in event.picks if "P" in pick.phase_hint]
            picksS = [pick for pick in event.picks if "S" in pick.phase_hint]
            vm = self.velocityModel
            source_z = depth*1e-3
            zGridMax = self.zGridMax
            receiversP = self.computeNewDistances(eventLatitude, eventLongitude, picksP, updatedStations)
            receiversS = self.computeNewDistances(eventLatitude, eventLongitude, picksS, updatedStations)
            node_intervals = self.node_intervals
            decimationFactor = self.decimationFactor
            extractedTTP = generateTT(vm, "P", source_z, zGridMax, receiversP, node_intervals, decimationFactor)
            extractedTTS = generateTT(vm, "S", source_z, zGridMax, receiversS, node_intervals, decimationFactor)
            for pick in event.picks:
                for pickP,newTT in zip(picksP, extractedTTP):
                    if pick.resource_id == pickP.resource_id:
                        pick.time = originTime + newTT
                for pickS,newTT in zip(picksS, extractedTTS):
                    if pick.resource_id == pickS.resource_id:
                        pick.time = originTime + newTT
        return catalog

    def loss(self, updatedStation, updatedCatalog):
        randomID = "".join(choices(al, k=6))
        self.writeCatalogFile(updatedCatalog, randomID)
        self.writeStationFile(updatedStation, randomID)
        self.writeVelocityFile(self.velocityModel, randomID)

    def runPSO(self):
        catalog = self.readEventFile()
        updatedStations = self.readStationFile(self.parDict["stationFile"])
        updatedCatalog = self.computeNewArrivals(catalog, updatedStations)
        self.writeHypoellipseConfigFile()
        self.loss(updatedStations, updatedCatalog)

    def plotResults(self):
        pass


# run application
myApp = main()
myApp.runPSO()