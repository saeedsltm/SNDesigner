from matplotlib.pyplot import title
import proplot as pplt
import numpy as np
import matplotlib.tri as tri
import cartopy.io.shapereader as shapereader
from shapely.geometry.polygon import Polygon
from matplotlib.path import Path

def extractIranBorder():
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
    triang = tri.Triangulation(x, y)
    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)
    maskedPoints = filterPointsInsideIran(Xi, Yi, polygon)
    zi = np.ma.array(zi, mask=maskedPoints)
    return zi

def filterPointsInsideIran(X, Y, poly):
    lons, lats = poly.exterior.coords.xy
    polygonPoints = [(x,y) for x,y in zip(lons, lats)]
    points = list(zip(X.flatten(), Y.flatten()))
    path = Path(polygonPoints, closed=True)
    inside = path.contains_points(points)
    return ~inside

def plotSeismicity(eventsX, eventsY, lonMin, lonMax, latMin, latMax, outName):
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
    
    axs[0].plot(eventsX, eventsY, marker="o", ms=5, mew=0.5, mfc="r", mec="k", ls="", zorder=10)
    fig.save("seismicity_{0:s}.png".format(outName))    

def plotGap(df, outName):
    df = df[(df.RMSs<2.0)&(df.ERHs>0)&(df.ERHs<10)&(df.ERZs>0)&(df.ERZs<10)]
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
    plotSeismicity(eventsX, eventsY, lonMin, lonMax, latMin, latMax, outName)
    for i,z in enumerate([eventsGap, eventsRMS, eventsERH, eventsERZ]):
        zi = interpolateData(eventsX, eventsY, xi, yi, z, iranPolygon)
        fig, axs = pplt.subplots(ncols=1, nrows=1, figwidth=5)
        # axs[0].pcolormesh(xi, yi ,zi,  cmap="greys")
        # axs[0].set_aspect("equal")
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
        axs[0].format(lonlim=(lonMin-1, lonMax+1), latlim=(latMin-1, latMax+1), labels=True)
        # axs[0].scatter(Xi.flatten(), Yi.flatten(), c=maskedPoints.astype("float").flatten(), cmap="jet", zorder=10)
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
        fig.colorbar(m, loc="b", label="{0:s} {1:s}".format(k, attributesDict[k][0]))
        fig.save("{0:s}_{1:s}.png".format(k, outName))