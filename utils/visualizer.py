import os
from copy import deepcopy

import cartopy.io.shapereader as shapereader
import matplotlib.tri as tri
import numpy as np
import proplot as pplt
from matplotlib.path import Path


def extractIranBorder():
    """Extract Iran boundray polygon

    Raises:
        ValueError: raise an error if "Iran" not found

    Returns:
        Polygon: a polygon contains Iran's border
    """
    countries = shapereader.natural_earth(
        resolution="110m",
        category="cultural",
        name="admin_0_countries")
    for country in shapereader.Reader(countries).records():
        if country.attributes["NAME"] == "Iran":
            iran = country.geometry
            return iran
    else:
        raise ValueError("Unable to find the Iran boundary.")


def interpolateData(x, y, xi, yi, z, polygon):
    """Interpolate given grid to new provided points and mask it by input polygon

    Args:
        x (1D-array): an array of points in x direction
        y (1D-array): an array of points in y direction
        xi (2D-array): 2D array of meshgrids in x direction
        yi (2D-array): 2D array of meshgrids in y direction
        z (1D-array): interpolated values
        polygon (Polygon): a polygon of Iran's border

    Returns:
        2-D array: 2-D array of interpolated points 
    """
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)
    maskedPoints = filterPointsInsideIran(Xi, Yi, polygon)
    zi = np.ma.array(zi, mask=maskedPoints)
    return zi


def filterPointsInsideIran(X, Y, poly):
    """Filter out data points if outside the given polygon

    Args:
        X (array): x points
        Y (array): y points
        poly (Polygon): desired polygon

    Returns:
        array.masked: an array contains masked points
    """
    lons, lats = poly.exterior.coords.xy
    polygonPoints = [(x, y) for x, y in zip(lons, lats)]
    points = list(zip(X.flatten(), Y.flatten()))
    path = Path(polygonPoints, closed=True)
    inside = path.contains_points(points)
    return ~inside


def plotSeismicity(eventsX, eventsY, lonMin, lonMax, latMin, latMax, outName):
    """Plot seismicity map

    Args:
        eventsX (array): events longitude
        eventsY (array): events latitude
        lonMin (float): min longitude of map
        lonMax (float): max longitude of map
        latMin (float): min latitude of map
        latMax (float): max latitude of map
        outName (str): map output name
    """
    fig, axs = pplt.subplots(ncols=1, nrows=1, figwidth=5)
    pplt.rc.reso = "med"
    fig, axs = pplt.subplots(nrows=1, refwidth=5, proj="cyl")
    axs.format(
        land=True,
        labels=True,
        coast=True,
        borders=True,
        landcolor="white",
        lonlines=5, latlines=5,
        facecolor="gray",
        xlabel="Longitude", ylabel="Latitude",
        suptitle="{n:d} events".format(n=eventsX.size)
    )
    axs[0].format(
        lonlim=(lonMin-1, lonMax+1),
        latlim=(latMin-1, latMax+1),
        labels=True
    )
    axs[0].plot(eventsX, eventsY, marker="o", ms=5,
                mew=0.5, mfc="r", mec="k", ls="", zorder=10)
    fig.save("seismicity_{0:s}.png".format(outName))


def plotStationDislocation(stationsDict, resultsPath, lonMin, lonMax, latMin, latMax):
    """Plot station dislocation

    Args:
        stationsDict (dict): a dictionary contains stations position
        resultsPath (str): results directory path
        lonMin (float): min longitude of map
        lonMax (float): max longitude of map
        latMin (float): min latitude of map
        latMax (float): max latitude of map
    """
    fig, axs = pplt.subplots(ncols=1, nrows=1, figwidth=5)
    pplt.rc.reso = "med"
    fig, axs = pplt.subplots(nrows=1, refwidth=5, proj="cyl")
    axs.format(
        land=True,
        labels=True,
        coast=True,
        borders=True,
        landcolor="white",
        lonlines=5, latlines=5,
        facecolor="gray",
        xlabel="Longitude", ylabel="Latitude",
        suptitle="station dislocation"
    )
    axs[0].format(
        lonlim=(lonMin-1, lonMax+1),
        latlim=(latMin-1, latMax+1),
        labels=True
    )
    initialStations = deepcopy(stationsDict)
    finalStations = np.loadtxt(os.path.join(resultsPath, "bestModel.dat"))
    for i, station in enumerate(list(initialStations.keys())):
        x = initialStations[station]["Lon"]
        y = initialStations[station]["Lat"]
        yy = finalStations[i]
        xx = finalStations[i+len(initialStations.keys())]
        axs[0].plot(x, y, "r^", ms=5, alpha=.5)
        axs[0].plot(xx, yy, "b^", ms=5, alpha=.5)
        axs[0].arrow(x, y, xx-x, yy-y, head_width=.05,
                     head_length=.05, zorder=10)
    fig.save("newStationLocations.png")


def plotResults(df, stationsDict, resultsPath, outName):
    """Plot final results

    Args:
        df (DataFrame): a dataframe contains relocated events 
        stationsDict (dict): a dictionary contains stations position
        resultsPath (str): results directory path
        outName (str): map output name
    """
    df = df[(df.RMSs < 2.0) & (df.ERHs > 0) & (
        df.ERHs < 10) & (df.ERZs > 0) & (df.ERZs < 10)]
    eventsX = df.LONs.values
    eventsY = df.LATs.values
    eventsGap = df.GAPs.values
    eventsRMS = df.RMSs.values
    eventsERH = df.ERHs.values
    eventsERZ = df.ERZs.values
    attributesDict = {
        # parameter    : [dimension, vmin, vmax]
        "Azimuthal gap": ["($\degree$)", 50, 360],
        "RMS": ["(s)", 0.0, 1.0],
        "Horizontal error": ["(km)", 0, 10],
        "Depth error": ["(km)", 0, 10]
    }
    N = 100
    iranPolygon = extractIranBorder()
    lonMin, latMin, lonMax, latMax = iranPolygon.bounds
    xi = np.linspace(lonMin-1, lonMax+1, N)
    yi = np.linspace(latMin-1, latMax+1, N)
    plotStationDislocation(stationsDict, resultsPath,
                           lonMin, lonMax, latMin, latMax)
    plotSeismicity(eventsX, eventsY, lonMin, lonMax, latMin, latMax, outName)
    for i, z in enumerate([eventsGap, eventsRMS, eventsERH, eventsERZ]):
        zi = interpolateData(eventsX, eventsY, xi, yi, z, iranPolygon)
        fig, axs = pplt.subplots(ncols=1, nrows=1, figwidth=5)
        pplt.rc.reso = "med"
        fig, axs = pplt.subplots(nrows=1, refwidth=5, proj="cyl")
        axs.format(
            land=True,
            labels=True,
            coast=True,
            borders=True,
            landcolor="white",
            lonlines=5, latlines=5,
            facecolor="gray",
        )
        k = list(attributesDict.keys())[i]
        axs[0].format(lonlim=(lonMin-1, lonMax+1),
                      latlim=(latMin-1, latMax+1), labels=True)
        m = axs[0].contourf(
            xi, yi, zi,
            cmap="solar_r",
            corner_mask=True,
            vmin=attributesDict[k][1],
            vmax=attributesDict[k][2],
            zorder=10)
        axs.format(
            xlabel="Longitude", ylabel="Latitude"
        )
        fig.colorbar(m, loc="b", label="{0:s} {1:s}".format(
            k, attributesDict[k][0]))
        fig.save("{0:s}_{1:s}.png".format(k, outName))
