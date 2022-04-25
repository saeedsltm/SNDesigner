import os
import warnings
from copy import deepcopy
from shutil import copy
from json import load
from multiprocessing import Pool
from pathlib import Path
from random import choices
from shutil import rmtree
from string import ascii_letters as al

import pykonal
from fstpso import FuzzyPSO
from LatLon import lat_lon as ll
from numpy import array, loadtxt, max, mean
from obspy import read_events
from obspy.geodetics.base import gps2dist_azimuth as gps

from utils.catalog2hypo import (catalog2hypoellipse,
                                createHypoellipseDefaultFile,
                                station2hypoellipse, velocity2hypoellipse)
from utils.extra import decoratortimer
from utils.hypo71 import runHypo71
from utils.hypoellipse import runHypoellipse
from  utils.hypocenter import runHypocenter
from utils.synthTime import extractTT, generateTTT
from utils.parseHypo import parseHypoellipseOutput, parseHypo71Output, parseHypocenterOutput
from utils.visualizer import plotGap


warnings.filterwarnings("ignore")


class main():
    def __init__(self):
        # reading main parameter file
        with open("parameters.json") as f:
            self.parDict = load(f)
        # reading "hypoellipse" defaults parameter file
        with open("hypoDefaults.json") as f:
            self.hypoDefaultsDict = load(f)
        # reading velocity file
        self.velocityModelDict = self.readVelocityFile(
            self.parDict["velocityFile"])
        # reading station file
        self.stationsDict = self.readStationFile(self.parDict["stationFile"])
        # setting "pykonal" parameters
        self.GridZMax = int(max(self.velocityModelDict["Z"]) + 1.0)
        self.node_intervals = self.parDict["velocityGridIntervals"]
        self.decimationFactor = self.parDict["decimationFactor"]
        # setting result directory
        self.resPath = os.path.join("results")
        # setting execution number
        self.run_id = 0
        # set global variable for storing traveltime tabels
        self.travelTimeDict = {}
        # reading catalog
        self.catalog = self.readEventFile()

    def makeResultDirecory(self):
        """make result directory
        """
        if os.path.exists(self.resPath):
            rmtree(self.resPath)
        self.resPath = Path(self.resPath)
        self.resPath.mkdir(parents=True, exist_ok=True)
        self.bestModelOut = os.path.join(self.resPath, "bestModel.dat")

    def readVelocityFile(self, stationFile):
        """reading velocity model from "STATION0.HYP" file 

        Args:
            stationFile (str): station file in "NORDIC" format

        Returns:
            (dict): a dictionary contains velocity model
        """
        emptyLines = 0
        velocityModelDict = {"Vp": [], "Z": [], "VpVs": 1.73, "Moho":46.0}
        with open(stationFile) as f:
            for l in f:
                if not l.strip():
                    emptyLines += 1
                if emptyLines == 2 and l.strip():
                    Vp, Z = [float(x) for x in l.split()[:2]]
                    velocityModelDict["Vp"].append(Vp)
                    velocityModelDict["Z"].append(Z)
                if emptyLines == 2 and len(l) > 20 and l[21] == "N":
                    _, Z = [float(x) for x in l.split()[:2]]
                    velocityModelDict["Moho"] = Z
                if emptyLines == 3 and l.strip():
                    VpVs = float(l[16:20])
                    velocityModelDict["VpVs"] = VpVs
                    break
        return velocityModelDict

    def readStationFile(self, stationFile):
        """read station information from "STATION0.HYP" file

        Args:
            stationFile (str): station file in "NORDIC" format

        Returns:
            dict: a dictionary contains stations information
        """
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
                    stations[code] = {"Lat": lat, "Lon": lon, "Elv": elv}
        return stations

    # @decoratortimer(2)
    def readEventFile(self):
        """read event file in "obspy" supported formats

        Returns:
            obspy.catalog: a catalog contains event information
        """
        catalog = read_events(self.parDict["eventsFile"])
        return catalog

    def generateTTTable(self):
        """generate travel time tables and store a bank of files
        """
        xMaxDist = self.parDict["xGridMax"] * self.parDict["velocityGridIntervals"][0]
        zMaxDist = self.parDict["zGridMax"] * self.parDict["velocityGridIntervals"][2]
        print("\n+++ Generating travel time tables for X range (0, {xMax:4.0f})km, and Z range (0, {zMax:4.0f})km\n".format(xMax=xMaxDist, zMax=zMaxDist))
        print("+++ Generating travel time tables for P phase ")
        generateTTT(
            self.velocityModelDict, "P",
            self.parDict["xGridMax"], self.parDict["zGridMax"],
            self.node_intervals, self.decimationFactor
        )
        print("+++ Generating travel time tables for S phase ")
        generateTTT(
            self.velocityModelDict, "S",
            self.parDict["xGridMax"], self.parDict["zGridMax"],
            self.node_intervals, self.decimationFactor
        )

    # @decoratortimer(2)
    def computeNewDistances(self, eventLat, eventLon, picks, stations):
        """given event coordinates it computes new station's distance

        Args:
            eventLat (float): event latitude in degree
            eventLon (float): event longitude in degree
            picks (obspy.picks): obspy picks for an event
            stations (dict): a dictionary contains stations information

        Returns:
            numpy.array: an array contains stations distance
        """
        stationsInPicks = [pick.waveform_id.station_code for pick in picks]
        stationLats = [stations[station]["Lat"] for station in stationsInPicks]
        stationLons = [stations[station]["Lon"] for station in stationsInPicks]
        distances = [gps(eventLat, eventLon, stationLat, stationLon)[
            0]*1e-3 for (stationLat, stationLon) in zip(stationLats, stationLons)]
        return [array([distance, 1.0, 0.0]) for distance in distances]

    # @decoratortimer(2)
    def writeCatalogFile(self, catalog, id):
        """write obspy catalog in hypoellipse format

        Args:
            catalog (obspy.catalog): an obspy catalog of events
            id (str): a unique id used for the file to be written 
        """
        catPath = Path(os.path.join("tmp"))
        catPath.mkdir(parents=True, exist_ok=True)
        outFile = os.path.join(catPath, "{id:6s}.pha".format(id=id))
        catalog2hypoellipse(catalog, outFile)

    # @decoratortimer(2)
    def writeStationFile(self, stations, id):
        """write station file in hypoellipse format

        Args:
            stations (list): a list contains stations latitude and longitude
            id (str): a unique id used for the file to be written

        Returns:
            dict: a dictionary contains stations information
        """
        stationsDict = deepcopy(self.stationsDict)
        for i, station in enumerate(stationsDict.keys()):
            stationsDict[station]["Lat"] = stations[i]
            stationsDict[station]["Lon"] = stations[i+int(len(stations)/2)]
        staPath = Path(os.path.join("tmp"))
        staPath.mkdir(parents=True, exist_ok=True)
        outFile = os.path.join(staPath, "{id:6s}.sta".format(id=id))
        station2hypoellipse(stationsDict, outFile)
        return stationsDict

    # @decoratortimer(2)
    def writeVelocityFile(self, velocity, id):
        """write velocity model file in hypoellipse format

        Args:
            velocity (dict): a dictionary contains velocity information
            id (str): a unique id used for the file to be written
        """
        velPath = Path(os.path.join("tmp"))
        velPath.mkdir(parents=True, exist_ok=True)
        outFile = os.path.join(velPath, "{id:6s}.prm".format(id=id))
        velocity2hypoellipse(velocity, outFile)

    def writeHypoellipseConfigFile(self):
        """write defaults parameter in "hypoellipse" format
        """
        defPath = Path(os.path.join("tmp"))
        defPath.mkdir(parents=True, exist_ok=True)
        outFile = os.path.join(defPath, "default.cfg")
        createHypoellipseDefaultFile(self.hypoDefaultsDict, outFile)

    # @decoratortimer(2)
    def computeNewArrivals(self, catalog, updatedStations):
        """given an obspy catalog and updated station locations it computes new arrival times

        Args:
            catalog (obspy.catalog): an obspy catalog contains events information
            updatedStations (dict): a dictionary contains stations information

        Returns:
            obspy.catalog: an updated obspy catalog
        """

        for event in catalog:
            preferred_origin = event.preferred_origin()
            originTime = preferred_origin.time
            eventLatitude = preferred_origin.latitude
            eventLongitude = preferred_origin.longitude
            depth = preferred_origin.depth
            picksP = [pick for pick in event.picks if "P" in pick.phase_hint]
            picksS = [pick for pick in event.picks if "S" in pick.phase_hint]
            source_z = depth*1e-3
            receiversP = self.computeNewDistances(
                eventLatitude, eventLongitude, picksP, updatedStations)
            receiversS = self.computeNewDistances(
                eventLatitude, eventLongitude, picksS, updatedStations)
            extractedTTP = []
            extractedTTS = []
            if len(receiversP):
                hdf5File = os.path.join(
                    "ttt", "dep{z:003d}{vt:s}.hdf5".format(z=int(source_z), vt="P"))
                if hdf5File not in self.travelTimeDict.keys():
                    traveltime = pykonal.fields.read_hdf(  # type: ignore
                        hdf5File)
                    self.travelTimeDict[hdf5File] = traveltime
                extractedTTP = extractTT(
                    self.travelTimeDict[hdf5File], receiversP)
            if len(receiversS):
                hdf5File = os.path.join(
                    "ttt", "dep{z:003d}{vt:s}.hdf5".format(z=int(source_z), vt="S"))
                if hdf5File not in self.travelTimeDict.keys():
                    traveltime = pykonal.fields.read_hdf(  # type: ignore
                        hdf5File)
                    self.travelTimeDict[hdf5File] = traveltime
                extractedTTS = extractTT(
                    self.travelTimeDict[hdf5File], receiversS)
            for pick in event.picks:
                for pickP, newTT in zip(picksP, extractedTTP):
                    if pick.resource_id == pickP.resource_id:
                        pick.time = originTime + newTT
                for pickS, newTT in zip(picksS, extractedTTS):
                    if pick.resource_id == pickS.resource_id:
                        pick.time = originTime + newTT
        return catalog

    def setSearchSpace(self, stationsDict):
        """setting bounds for stations coordinates

        Args:
            stationsDict (dict): a dictionary contains stations information

        Returns:
            list: a list contains each station min-max bounds
        """
        stationLats = [stationsDict[station]["Lat"]
                       for station in stationsDict]
        stationLons = [stationsDict[station]["Lon"]
                       for station in stationsDict]
        incLatMin = self.parDict["stationLatBound"][0]
        incLatMax = self.parDict["stationLatBound"][1]
        incLonMin = self.parDict["stationLatBound"][0]
        incLonMax = self.parDict["stationLatBound"][1]
        boundsLat = [[lat + incLatMin, lat + incLatMax] for lat in stationLats]
        boundsLon = [[lon + incLonMin, lon + incLonMax] for lon in stationLons]
        searchSpace = boundsLat + boundsLon
        return searchSpace

    def writeFSTPSOResults(self, bestStations):
        """write fst-pso results in a file

        Args:
            bestStations (list): a list contains best stations coordinates
        """
        with open(self.bestModelOut, "a") as f:
            f.write(" ".join("{0:7.2f}".format(e) for e in bestStations[0].X))
            f.write("\n")
        stationsDict = deepcopy(self.stationsDict)
        bestStations = [v for v in bestStations[0].X]
        for i, station in enumerate(stationsDict.keys()):
            stationsDict[station]["Lat"] = bestStations[i]
            stationsDict[station]["Lon"] = bestStations[i +
                                                        int(len(bestStations)/2)]
        outFile = os.path.join(self.resPath, "final.sta")
        station2hypoellipse(stationsDict, outFile)

    def evaluateMisfit(self, metric, gap, rms, erh, erz):
        """evalute misfit value based on location accuracy

        Args:
            metric (str): metric name ("GAP", "RMS", "ERH", "ERZ")
            gap (list): a list contains all events azimuthal gap (degree) 
            rms (list): a list contains all events rms (s) 
            erh (list): a list contains all events horizontal error (km) 
            erz (list): a list contains all events depth error (km) 

        Returns:
            float: a misfit value
        """
        misfits = {
            "GAP": mean(gap),
            "RMS": mean(rms),
            "ERH": mean(erh),
            "ERZ": mean(erz)
        }
        return misfits[metric]

    # @decoratortimer(2)
    def loss(self, updatedStations):
        """given an updated location of stations it computes misfit based on event location accuracy

        Args:
            updatedStations (dict): a dictionary contains stations information

        Returns:
            float: loss value
        """
        randomID = "".join(choices(al, k=6))
        stationsDict = self.writeStationFile(updatedStations, randomID)
        updatedCatalog = self.computeNewArrivals(self.catalog, stationsDict)
        self.writeCatalogFile(updatedCatalog, randomID)
        self.writeVelocityFile(self.velocityModelDict, randomID)
        root = os.getcwd()
        os.chdir("tmp")
        if self.parDict["misfitFunction"] == "hypoellipse":
            GAP, RMS, ERH, ERZ = runHypoellipse(randomID)
            misfit = self.evaluateMisfit(
                self.parDict["metric"], GAP, RMS, ERH, ERZ)
            cmd = "rm {prefix}*".format(prefix=randomID)
            os.system(cmd)
            os.chdir(root)
            return misfit
        elif self.parDict["misfitFunction"] == "hypo71":
            GAP, RMS, ERH, ERZ = runHypo71(randomID, self.hypoDefaultsDict)
            misfit = self.evaluateMisfit(
                self.parDict["metric"], GAP, RMS, ERH, ERZ)
            cmd = "rm {prefix}*".format(prefix=randomID)
            os.system(cmd)
            os.chdir(root)
            return misfit

    def lossparallel(self, particles):
        """execute loos function in parallel

        Args:
            particles (list): a nested list contains new stations location

        Returns:
            list: a list of loss values
        """
        if self.parDict["multiProcessingMode"] == 0:
            numberOfCPU = 1
        else:
            numberOfCPU = os.cpu_count()
        p = Pool(numberOfCPU)
        results = p.map(self.loss, particles)
        p.close()
        return results

    def runPSO(self):
        """main function of fst-pso
        """
        self.makeResultDirecory()
        self.writeHypoellipseConfigFile()
        FP = FuzzyPSO(logfile=os.path.join(self.resPath, "log.dat"))
        searchSpace = self.setSearchSpace(self.stationsDict)
        FP.set_search_space(searchSpace)
        bestFitnessOutFile = os.path.join(
            self.resPath, "bestFitness_{0:d}.dat".format(self.run_id+1))
        if self.parDict["numberOfModel"] != 0:
            FP.set_swarm_size(self.parDict["numberOfModel"])
        if self.parDict["multiProcessingMode"] != 0:
            FP.set_parallel_fitness(self.lossparallel)
            if self.parDict["numberOfIteration"] != 0:
                result = FP.solve_with_fstpso(
                    max_iter=self.parDict["numberOfIteration"],
                    max_iter_without_new_global_best=self.parDict["maxIterWithoutNewGlobalBest"],
                    dump_best_fitness=bestFitnessOutFile)
            else:
                result = FP.solve_with_fstpso(
                    max_iter_without_new_global_best=self.parDict["maxIterWithoutNewGlobalBest"],
                    dump_best_fitness=bestFitnessOutFile)
        else:
            FP.set_fitness(self.loss, skip_test=True)
            if self.parDict["numberOfIteration"] != 0:
                result = FP.solve_with_fstpso(
                    max_iter=self.parDict["numberOfIteration"],
                    max_iter_without_new_global_best=self.parDict["maxIterWithoutNewGlobalBest"],
                    dump_best_fitness=bestFitnessOutFile)
            else:
                result = FP.solve_with_fstpso(
                    max_iter_without_new_global_best=self.parDict["maxIterWithoutNewGlobalBest"],
                    dump_best_fitness=bestFitnessOutFile)
        self.writeFSTPSOResults(result)

    def finalLocation(self, hypoellipse=False, hypo71=False, hypocenter=True):
        outName = "finRun"
        self.writeCatalogFile(self.catalog, outName)
        self.writeVelocityFile(self.velocityModelDict, outName)
        copy(os.path.join("results", "final.sta"), os.path.join("tmp", "{0:s}.sta".format(outName)))
        if hypoellipse:
            root = os.getcwd()
            os.chdir("tmp")
            runHypoellipse(outName)
            os.chdir(root)
            return os.path.join("tmp", "{0:s}.out".format(outName))
        elif hypo71:
            root = os.getcwd()
            os.chdir("tmp")            
            runHypo71(outName, self.hypoDefaultsDict)
            os.chdir(root)
            return os.path.join("tmp", "{0:s}_h71.out".format(outName))
        elif hypocenter:
            root = os.getcwd()
            os.chdir("tmp")            
            runHypocenter(outName, self.velocityModelDict, self.stationsDict, self.hypoDefaultsDict)
            os.chdir(root)
            return os.path.join("tmp", "{0:s}.out".format(outName))

    @decoratortimer(2)
    def plotResults(self):
        """simple plot of the results
        """
        import matplotlib.pyplot as plt
        initialStations = deepcopy(self.stationsDict)
        finalStations = loadtxt(os.path.join(self.resPath, "bestModel.dat"))
        for i, station in enumerate(list(initialStations.keys())):
            x = initialStations[station]["Lon"]
            y = initialStations[station]["Lat"]
            yy = finalStations[i]
            xx = finalStations[i+len(initialStations.keys())]
            plt.plot(x, y, "r^", ms=5, alpha=.5)
            plt.plot(xx, yy, "b^", ms=5, alpha=.5)
            plt.arrow(x, y, xx-x, yy-y, head_width=.05,
                      head_length=.05, zorder=10)
        plt.savefig(os.path.join(self.resPath, "BestModel.png"))


# run application
myApp = main()
# myApp.generateTTTable()
# myApp.runPSO()
# myApp.plotResults()
dfHypocenter = parseHypocenterOutput("events/select.out")
plotGap(dfHypocenter, "initial")
hypo71Out = myApp.finalLocation(hypo71=True)
dfHypo71 = parseHypo71Output(hypo71Out)
plotGap(dfHypo71, "finHypo71")
