import proplot as pplt
import numpy as np
import matplotlib.tri as tri
import cartopy.io.shapereader as shapereader
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from matplotlib.path import Path
import matplotlib.pyplot as plt

def extractIranBorder():
##    iran = Polygon()
    countries = shapereader.natural_earth(
        resolution="110m",
        category="cultural",
        name="admin_0_countries")
    # Find the Iran boundary polygon
    for country in shapereader.Reader(countries).records():
        if country.attributes["NAME"] == "Iran":
            iran = country.geometry
            return iran
    else:
        raise ValueError("Unable to find the Iran boundary.")


def filterPointsInsideIran(X, Y, poly):
    lons, lats = poly.exterior.coords.xy
    polygonPoints = [(x,y) for x,y in zip(lons, lats)]
    points = list(zip(X.flatten(), Y.flatten()))
    path = Path(polygonPoints, closed=True)
    inside = path.contains_points(points)
    return inside.reshape(X.shape)


LatMin = 24.00
LatMax = 42.00
LonMin = 42.00
LonMax = 64.00
N = 10
xi = np.linspace(LonMin, LonMax, N)
yi = np.linspace(LatMin, LatMax, N)
Xi, Yi = np.meshgrid(xi, yi)
inside = filterPointsInsideIran(Xi, Yi, extractIranBorder())

plt.scatter(Xi.flatten(), Yi.flatten(), c=inside.astype("float"))
plt.show()
