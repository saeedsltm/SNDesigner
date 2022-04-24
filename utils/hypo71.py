import os
import platform

from LatLon import lat_lon as ll
from utils.extra import mean


def addHypo71Defaults(fileObject):
    """Add RESET TESTS to hypo71 input file

    Args:
        fileObject (fileObject): file object of hypo71 input file
    """
    fileObject.write(
        "HEAD                     GENERATED USING 'hypo71.py' CODE\n")
    fileObject.write("RESET TEST(01)=0.1\n")
    fileObject.write("RESET TEST(02)=10.\n")
    fileObject.write("RESET TEST(03)=1.5\n")
    fileObject.write("RESET TEST(04)=0.05\n")
    fileObject.write("RESET TEST(05)=5.\n")
    fileObject.write("RESET TEST(06)=4.\n")
    fileObject.write("RESET TEST(11)=40.\n")
    fileObject.write("\n")


def addHypo71Stations(fileObject, stationDict):
    """Add stations to hypo71 input file

    Args:
        fileObject (fileObject): file object of hypo71 input file
        stationDict (dict): a dictionary contains stations information
    """
    for station in stationDict:
        fileObject.write(stationDict[station])
    fileObject.write("\n")


def addHypo71VelocityModel(fileObject, velocityModelDict):
    """Add velocity model to hypo71 input file

    Args:
        fileObject (fileObject): file object of hypo71 input file
        velocityModelDict (dict): a dictionary contains P velocity and Vp/Vs ratio
    """
    Vp, Z, _ = velocityModelDict.items()
    for v, z in zip(Vp[1], Z[1]):
        layer = "  {0:5.3f} {1:6.3f}\n".format(v, z)
        fileObject.write(layer)
    fileObject.write("\n")


def addHypo71ControlLine(fileObject, controlLine):
    """Add control line in hypo71 input file

    Args:
        fileObject (fileObject): file object of hypo71 input file
        controlLine (str): a string of hypo71 control line
    """
    fileObject.write(controlLine)


def addHypo71Phase(objectFile, phaseFile):
    """Add phases to hypo71 input file

    Args:
        objectFile (fileObject): file object of hypo71 input file
        phaseFile (list): a list contains hypo71 phase line
    """
    with open(phaseFile) as f:
        for l in f:
            objectFile.write(l)


def parseHypoellipseStation(hypoellipseStationFile):
    """Parse hypoellipse station file

    Args:
        hypoellipseStationFile (str): station file in hypoellipse format

    Returns:
        dict: a dictionary contains stations information
    """
    stationDict = {}
    with open(hypoellipseStationFile) as f:
        for l in f:
            if "*" not in l:
                staCode = l[:4].strip()
                latDeg, latMin = int(l[4:6]), float(l[7:12])
                lat = ll.Latitude(latDeg, latMin)
                lonDeg, lonMin = int(l[13:16]), float(l[17:22])
                lon = ll.Longitude(lonDeg, lonMin)
                elv = int(l[22:27])
                stationDict[staCode] = "  {0:4s}{1:2d}{2:05.2f}N {3:2d}{4:05.2f}E{5:4d}\n".format(
                    staCode, int(lat.decimal_degree), lat.decimal_minute, int(lon.decimal_degree), lon.decimal_minute, elv)
    return stationDict


def parseHypoellipseVelocity(hypoellipseVelocityModelFile):
    """Parse hypoellipse velocity file

    Args:
        hypoellipseVelocityModelFile (str): velocity file in hypoellipse format

    Returns:
        dict: a dictionary contains velocity model
    """
    velocityModelDict = {"Vp": [], "Z": [], "VpVs": []}
    with open(hypoellipseVelocityModelFile) as f:
        for l in f:
            velocityModelDict["Vp"].append(float(l.split()[1]))
            velocityModelDict["Z"].append(float(l.split()[2]))
            velocityModelDict["VpVs"].append(float(l.split()[3]))
    return velocityModelDict


def prepareHypo71ControlLine(hypoellipseVelocityModelFile, defaultsDict):
    """_summary_

    Args:
        hypoellipseVelocityModelFile (str): velocity model fil in hypoellipse format
        defaultsDict (dict): a dictionary contains defaults for hypo71 control line

    Returns:
        str: hypo71 control line
    """  
    velocityModelDict = parseHypoellipseVelocity(hypoellipseVelocityModelFile)
    VpVs = mean(velocityModelDict["VpVs"])
    controlLine = "{0:4.0f}.{1:4.0f}.{2:4.0f}.{3:5.2f}    4    0    0    1    1    0    0 0111\n".format(
        defaultsDict["startingDepth"], defaultsDict["distanceWeighting"][1], defaultsDict["distanceWeighting"][2], VpVs)
    return controlLine


def getStatistic(rootName):
    """Calculate statistic results

    Args:
        rootName (str): root name of hypo71 output files

    Returns:
        tuple: a tuple contains lists of Azimuthal gap, RMS, Horizontal and depth errors
    """
    GAP, RMS, ERH, ERZ = [], [], [], []
    with open("{0}_h71.out".format(rootName)) as f:
        next(f)
        for l in f:
            if "*" not in l:
                try:
                    GAP.append(float(l[54:57]))
                    RMS.append(float(l[62:67]))
                    ERH.append(float(l[67:72]))
                    ERZ.append(float(l[72:77]))
                except ValueError:
                    continue
    return GAP, RMS, ERH, ERZ


def createHypo71PhaseFile(rootName, defaultsDict):
    """Create phase file in hypo71 format

    Args:
        rootName (str): root name of hypo71 output files
    """
    hypo71PhaseFile = "{0}_h71.pha".format(rootName)
    hypoellipseStationFile = "{0}.sta".format(rootName)
    hypoellipseVelocityModelFile = "{0}.prm".format(rootName)
    hypoellipsePhaseFile = "{0}.pha".format(rootName)
    with open(hypo71PhaseFile, "w") as f:
        addHypo71Defaults(f)
        stationDict = parseHypoellipseStation(hypoellipseStationFile)
        addHypo71Stations(f, stationDict)
        velocityModel = parseHypoellipseVelocity(hypoellipseVelocityModelFile)
        addHypo71VelocityModel(f, velocityModel)
        controlLine = prepareHypo71ControlLine(hypoellipseVelocityModelFile, defaultsDict)
        addHypo71ControlLine(f, controlLine)
        addHypo71Phase(f, hypoellipsePhaseFile)


def createHypo71InputFile(rootName):
    """Create input file in hypo71 format

    Args:
        rootName (str): root name of hypo71 output files
    """
    outName = "{0}.inp".format(rootName)
    with open(outName, "w") as f:
        f.write("{0}_h71.pha\n".format(rootName))
        f.write("{0}_h71.prt\n".format(rootName))
        f.write("{0}_h71.out\n".format(rootName))
        f.write("{0}_h71.res\n".format(rootName))
        f.write("\n\n")


def runHypo71(rootName, defaultsDict):
    """Run hypo71 program to locate earthquakes

    Args:
        rootName (str): root name of hypo71 output files
        defaultsDict (dict): a dictionary contains default values for hypo71 TESTs

    Returns:
        tuple: a tuple contains lists of Azimuthal gap, RMS, Horizontal and depth errors
    """
    createHypo71PhaseFile(rootName, defaultsDict)
    createHypo71InputFile(rootName)
    if platform.system() == "Linux":
        cmd = "../utils/hypo71Main < {0}.inp > /dev/null".format(rootName)
        os.system(cmd)
    elif platform.system() == "Windows":
        cmd = "..\\utils\\hypo71Main.exe < {0}.inp > NUL".format(rootName)
        os.system(cmd)
    GAP, RMS, ERH, ERZ = getStatistic(rootName)
    for i in ["{0}_h71.pha".format(rootName),
              "{0}_h71.out".format(rootName),
              "{0}_h71.prt".format(rootName),
              "{0}_h71.res".format(rootName),
              "{0}_h71.prm".format(rootName),
              "{0}_h71.sta".format(rootName),
              "{0}.inp".format(rootName)]:
        if os.path.exists(i):
            os.remove(i)
    return GAP, RMS, ERH, ERZ
