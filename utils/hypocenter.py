import os
import platform
import sys
from LatLon import lat_lon as ll
from shutil import copy

from utils.extra import decoratortimer, mean

# @decoratortimer(2)
def prepareStationFile(velocityModelDict, stationsDict, defaultsDict):
    with open(os.path.join("..", "events", "STATION0.HYP")) as f, open("STATION0.HYP", "w") as g:
        for l in f:
            if l.startswith("RESET"):
                g.write(l)
        g.write("\n")
        for station in stationsDict.keys():
            lat = ll.Latitude(stationsDict[station]["Lat"])
            lon = ll.Longitude(stationsDict[station]["Lon"])
            elv = stationsDict[station]["Elv"]
            stationLine = "  {code:4s}{latDeg:2.0f}{latMin:05.2f}{latHem:1s} {lonDeg:2.0f}{lonMin:05.2f}{lonHem:1s}{elv:00004.0f}\n".format(
                code=station,
                latDeg=lat.degree, latMin=lat.decimal_minute, latHem=lat.get_hemisphere(),
                lonDeg=lon.degree, lonMin=lon.decimal_minute, lonHem=lon.get_hemisphere(),
                elv=elv
            )
            g.write(stationLine)
        g.write("\n")
        for v, z in zip(velocityModelDict["Vp"], velocityModelDict["Z"]):
            velocityLayer = "  {v:4.2f}     {z:4.1f}            \n".format(v=v, z=z)
            if z == velocityModelDict["Moho"]:
                velocityLayer = "  {v:4.2f}     {z:4.1f}      N     \n".format(v=v, z=z)
            g.write(velocityLayer)
        g.write("\n")
        controlLine = "{startingDepth:4.1f} {xNear:4.0f}.{xFar:4.0f}. {vpvs:4.2f}    \n".format(
            startingDepth=defaultsDict["startingDepth"],
            xNear=defaultsDict["distanceWeighting"][1],
            xFar=defaultsDict["distanceWeighting"][2],
            vpvs=velocityModelDict["VpVs"]
        )
        g.write(controlLine)
        g.write("BIN")

def preparePhaseFile():
    copy(os.path.join("..", "events", "select.out"), "select.out")

def runHypocenter(rootName, velocityModelDict, stationsDict, defaultsDict):
    preparePhaseFile()
    prepareStationFile(velocityModelDict, stationsDict, defaultsDict)
    with open("hypocenter.inp", "w") as f:
        f.write("select.out\n")
        f.write("n\n")
    cmd = "hyp < hypocenter.inp > /dev/null"
    os.system(cmd)
    copy("hyp.out", "{0:s}.out".format(rootName))
    return os.path.join("{0:s}.out".format(rootName))
