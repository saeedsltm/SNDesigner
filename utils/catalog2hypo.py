from LatLon import lat_lon as ll


def upper(s):
    """make the input string uppercase

    Args:
        s (str): input string

    Returns:
        str: output string
    """
    return s.upper()


def catalog2hypoellipse(catalog, outFile):
    """convert an obspy catalog to "hypoellipse" format

    Args:
        catalog (obspy.catalog): an obspy catalog
        outFile (str): name of output file which will be saved in
    """
    with open(outFile, "w") as f:
        for event in catalog:
            picks = {}
            for pick in event.picks:
                station = pick.waveform_id.station_code
                if station not in picks:
                    picks[station] = {"P": {"phase": " ", "time": "               ", "w": " "}, "S": {
                        "phase": " ", "time": "      ", "w": " "}}
                if "P" in pick.phase_hint:
                    try:
                        weight = pick.extra.get("nordic_pick_weight")["value"]
                    except AttributeError:
                        weight = " "
                    picks[station]["P"]["phase"] = pick.phase_hint
                    picks[station]["P"]["time"] = pick.time
                    picks[station]["P"]["w"] = weight
                elif "S" in pick.phase_hint:
                    try:
                        weight = pick.extra.get("nordic_pick_weight")["value"]
                    except AttributeError:
                        weight = " "
                    picks[station]["S"]["phase"] = pick.phase_hint
                    picks[station]["S"]["time"] = pick.time
                    picks[station]["S"]["w"] = weight
            for station in picks:
                if not picks[station]["P"]["phase"][0].strip():
                    continue
                phaseP = picks[station]["P"]["phase"][0].upper()
                arrivaltimeP = picks[station]["P"]["time"].strftime("%y%m%d%H%M%S.%f")[
                    :15]
                wP = picks[station]["P"]["w"]
                if picks[station]["S"]["phase"].strip():
                    phaseS = picks[station]["S"]["phase"][0].upper()
                    arrivaltimeS = picks[station]["S"]["time"] - \
                        picks[station]["P"]["time"].replace(
                            second=0, microsecond=0)
                    arrivaltimeS = "{arrivaltimeS:6.2f}".format(
                        arrivaltimeS=arrivaltimeS)
                    wS = picks[station]["S"]["w"]
                else:
                    phaseS = picks[station]["S"]["phase"]
                    arrivaltimeS = picks[station]["S"]["time"]
                    wS = picks[station]["S"]["w"]
                f.write("{station:4s} {phaseP:1s} {wP:1s} {arrivaltimeP:15s}      {arrivaltimeS:6s} {phaseS:1s} {wS:1s}          \n".format(
                    station=station, phaseP=phaseP, arrivaltimeP=arrivaltimeP, wP=wP, phaseS=phaseS, arrivaltimeS=arrivaltimeS, wS=wS))
            f.write("                 10\n")


def station2hypoellipse(stationDict, outFile):
    """convert stations to "hypoellipse" format

    Args:
        stationDict (dict): a dictionary contains stations information
        outFile (str): name of output file which will be saved in
    """
    with open(outFile, "w") as f:
        for sta in sorted(stationDict.keys(), key=lambda x: (len(x), upper(x), sorted(x))):
            if "*" in sta:
                continue
            lat = ll.Latitude(stationDict[sta]["Lat"])
            lon = ll.Longitude(stationDict[sta]["Lon"])
            elv = stationDict[sta]["Elv"]
            f.write("{station:4s}{latDeg:2.0f}{latHemispher:1s}{latDec:5.2f} {lonDeg:3.0f}{lonHemispher:1s}{lonDec:5.2f} {elevation:4.0f}\n".format(
                station=sta, latDeg=lat.degree, latHemispher=lat.get_hemisphere(), latDec=lat.decimal_minute, lonDeg=lon.degree, lonHemispher=lon.get_hemisphere(), lonDec=lon.decimal_minute, elevation=elv
            ))
            f.write("{station:4s}*     0     1.00\n".format(station=sta))


def velocity2hypoellipse(velocityDict, outFile):
    """Convert velocity to "hypoellipse" format

    Args:
        velocityDict (dict): a dictionary contains P velocity and Vp/Vs ratio
        outFile (str): name of output file which will be saved in
    """
    with open(outFile, "w") as f:
        for v, z in zip(velocityDict["Vp"], velocityDict["Z"]):
            f.write("VELOCITY             {Vp:4.2f} {Z:5.2f} {VpVs:4.2f}\n".format(
                Vp=v, Z=z, VpVs=velocityDict["VpVs"]))


def createHypoellipseDefaultFile(defaultsDict, outFile):
    """Creating hypoellipse defaults file

    Args:
        defaultsDict (dict): a dictionary contains hypoellipse defaults
        outFile (str): name of output file which will be saved in
    """
    with open(outFile, "w") as f:
        f.write("reset test         1    {0:6.2f}\n".format(
            defaultsDict["VpVs"]))
        f.write("reset test         2    {0:6.2f}\n".format(
            defaultsDict["elevationCorrection"][0]))
        f.write("reset test         8    {0:6.2f}\n".format(
            defaultsDict["elevationCorrection"][1]))
        f.write("reset test         5    {0:6.2f}\n".format(
            defaultsDict["startingDepth"]))
        f.write("reset test         10   {0:6.2f}\n".format(
            defaultsDict["distanceWeighting"][0]))
        f.write("reset test         11   {0:6.2f}\n".format(
            defaultsDict["distanceWeighting"][1]))
        f.write("reset test         12   {0:6.2f}\n".format(
            defaultsDict["distanceWeighting"][2]))
        f.write("reset test         21   {0:6.2f}\n".format(
            defaultsDict["maximumNumberOfIterations"]))
        f.write("reset test         29   {0:6.2f}\n".format(
            defaultsDict["standardErrorForArrivalTimesWithWeightCode0"]))
        f.write("reset test         38   {0:6.2f}\n".format(
            defaultsDict["locateWithS"]))
        f.write("reset test         39   {0:6.2f}\n".format(
            defaultsDict["factorForWeightsOfSAndS_PTimes"]))
        f.write("summary option     {0:1.0f}\n".format(
            defaultsDict["summaryOption"]))
        f.write("printer option     {0:1.0f}\n".format(
            defaultsDict["printerOption"]))
        f.write("constants noprint  {0:1.0f}\n".format(
            defaultsDict["constantsNoPrint"]))
        f.write("compress option    {0:1.0f}\n".format(
            defaultsDict["compressOption"]))
        f.write("tabulation option  {0:1.0f}\n".format(
            defaultsDict["tabulationOption"]))
        f.write("weight option      {0:4.2f} {0:4.2f} {0:4.2f}\n".format(
            defaultsDict["weightOption"][0], defaultsDict["weightOption"][1], defaultsDict["weightOption"][2]))
        f.write("ignore summary rec {0:1.0f}\n".format(
            defaultsDict["ignoreSummaryRec"]))
        f.write("header option      {0:s}\n".format(
            defaultsDict["headerOption"]))
