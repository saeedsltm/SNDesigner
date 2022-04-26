from copy import deepcopy
from datetime import datetime as dt
from math import sqrt

import pandas as pn
from LatLon import lat_lon as ll
from obspy import read_events
from obspy.geodetics import degrees2kilometers as d2k

headers = {"OTs": [],
           "LATs": [],
           "LONs": [],
           "DEPs": [],
           "GAPs": [],
           "RMSs": [],
           "ERHs": [],
           "ERZs": []
           }


def parseOT(l, hypo71=False):
    """Parse origin time line

    Args:
        l (str): a line string
        hypo71 (bool, optional): choose input line format. Defaults to False.

    Returns:
        datetime: an origin time
    """
    dateTimeIndex = [5, 7, 10, 11, 12, 13, 15, 16, 18, 19]
    fmt = " %Y%m%d %H%M %S.%f"
    reclen = 20
    if hypo71:
        dateTimeIndex = [0, 1, 2, 3, 4, 7, 8, 9, 10, 12, 13, 15, 16]
        fmt = "%y%m%d %H%M %S.%f"
        reclen = 17
    for i in dateTimeIndex:
        l = l[:i] + l[i].replace(" ", "0") + l[i+1:]
    ot = dt.strptime(l[:reclen], fmt)
    return ot


def parseStatistics(l):
    """Parse earthquake information and statistics

    Args:
        l (str): a line string

    Returns:
        tuple: a tuple contains earthquake information
    """
    lat = ll.Latitude(degree=float(l[21:23]),
                      minute=float(l[24:29])).decimal_degree
    lon = ll.Longitude(degree=float(l[31:33]),
                       minute=float(l[34:39])).decimal_degree
    dep = float(l[41:46])
    gap = float(l[60:63])
    rms = float(l[66:70])
    return lat, lon, dep, gap, rms


def parseOriginLine(l):
    """Parse earthquake information

    Args:
        l (str): a line string

    Returns:
        tuple: a tuple contains earthquake information
    """
    ot = parseOT(l)
    lat, lon, depth, gap, rms = parseStatistics(l)
    return ot, lat, lon, depth, gap, rms


def parseUncertaintyLine(l):
    """Parse earthquake uncertainty

    Args:
        l (str): a line string

    Returns:
        tuple: a tuple contains earthquake uncertainty
    """
    seh = float(l[3:7])
    sez = float(l[8:12])
    return seh, sez


def parseHypoellipseOutput(hypOut):
    """Parse hypoellipse output

    Args:
        hypOut (str): path to hypoellipse output file

    Returns:
        DataFrame: a pandas DataFrame include all events information
    """
    df = pn.DataFrame(headers)
    with open(hypOut) as f:
        ot, lat, lon, dep, gap, rms, seh, sez = [None]*8
        for l in f:
            if "date    origin" in l:
                ot, lat, lon, dep, gap, rms = parseOriginLine(next(f))
            elif "seh  sez q sqd" in l:
                seh, sez = parseUncertaintyLine(next(f))
                eventDict = deepcopy(headers)
                keys = eventDict.keys()
                values = [ot, lat, lon, dep, gap, rms, seh, sez]
                for k, v in zip(keys, values):
                    eventDict[k].append(v)
                df = pn.concat([df, pn.DataFrame(eventDict)])
    return df


def parseHypo71Output(hypOut):
    """Parse hypo71 output

    Args:
        hypOut (str): path to hypoellipse output file

    Returns:
        DataFrame: a pandas DataFrame include all events information
    """
    df = pn.DataFrame(headers)
    ot, lat, lon, dep, gap, rms, seh, sez = [None]*8
    with open(hypOut) as f:
        next(f)
        for l in f:
            try:
                ot = parseOT(l, hypo71=True)
                lat = ll.Latitude(degree=float(
                    l[18:20]), minute=float(l[21:26])).decimal_degree
                lon = ll.Longitude(degree=float(
                    l[28:30]), minute=float(l[31:36])).decimal_degree
                dep = float(l[38:43])
                gap = float(l[54:57])
                rms = float(l[62:67])
                seh = float(l[67:72])
                sez = float(l[72:77])
                eventDict = deepcopy(headers)
                keys = eventDict.keys()
                values = [ot, lat, lon, dep, gap, rms, seh, sez]
                for k, v in zip(keys, values):
                    eventDict[k].append(v)
                df = pn.concat([df, pn.DataFrame(eventDict)])
            except ValueError:
                continue
    return df


def parseHypocenterOutput(hypOut):
    """Parse hypocenter output

    Args:
        hypOut (str): path to hypoellipse output file

    Returns:
        DataFrame: a pandas DataFrame include all events information
    """
    catalog = read_events(hypOut)
    df = pn.DataFrame(headers)
    ot, lat, lon, dep, gap, rms, seh, sez = [None]*8
    for event in catalog:
        preferred_origin = event.preferred_origin()
        if not preferred_origin.quality:
            continue
        eventDict = deepcopy(headers)
        ot = preferred_origin.time
        lat = preferred_origin.latitude
        lon = preferred_origin.longitude
        dep = preferred_origin.depth*1e-3
        rms = preferred_origin.quality.standard_error
        gap = preferred_origin.quality.azimuthal_gap
        seh = d2k(sqrt(preferred_origin.latitude_errors.uncertainty **
                  2 + preferred_origin.longitude_errors.uncertainty**2))
        sez = preferred_origin.depth_errors.uncertainty*1e-3
        keys = eventDict.keys()
        values = [ot, lat, lon, dep, gap, rms, seh, sez]
        for k, v in zip(keys, values):
            eventDict[k].append(v)
        df = pn.concat([df, pn.DataFrame(eventDict)])
    return df
