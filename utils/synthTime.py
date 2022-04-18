import os
from pathlib import Path

import pykonal
from numpy import array, float16, ones
from tqdm.contrib import tzip

from utils.extra import decoratortimer


# @decoratortimer(2)
def generateVelocityGrid(velocityModelDict, nXnYnZ, node_intervals=(1., 1., 1.), min_coords=(0., 0., 0.), df=0):
    """Generating velocity grid for both P and S velocities

    Args:
        velocityModelDict (dict): a dictionary contains P velocity and Vp/Vs ration
        nXnYnZ (tuple): a tuple contains number of nodes in X,Y and Z directions
        node_intervals (tuple, optional): node intervals in X,Y and Z directions. Defaults to (1., 1., 1.).
        min_coords (tuple, optional): starting point of velocity grid. Defaults to (0., 0., 0.).
        df (int, optional): decimation factor. Defaults to 0.

    Raises:
        RuntimeError: when velocity interface (?km) exceeds number of grids (?km) in Z direction!

    Returns:
        tuple: computed velocity grids for P and S
    """
    DATA_TYPE = float16
    vgP = ones((nXnYnZ), dtype=DATA_TYPE)
    maxZ = nXnYnZ[2] * node_intervals[2]
    Vp, Z, VpVs = [v for _, v in velocityModelDict.items()]
    for vp, z in zip(Vp, Z):
        indx = int(z * 1/node_intervals[2])
        if z >= maxZ:
            errMsg = "Velocity interface ({0}km) exceeds number of grids ({1}km) in Z direction!".format(
                z, maxZ)
            raise RuntimeError(errMsg)
        vgP[:, :, indx:] = vp
    df = 2**df
    VG = vgP[::df, :, ::df]
    node_intervals = array(node_intervals) * df
    # P velocity
    solverP = pykonal.EikonalSolver(coord_sys="cartesian")
    solverP.velocity.min_coords = min_coords
    solverP.velocity.node_intervals = node_intervals
    solverP.velocity.npts = VG.shape
    solverP.velocity.values = VG
    # S velocity
    solverS = pykonal.EikonalSolver(coord_sys="cartesian")
    solverS.velocity.min_coords = min_coords
    solverS.velocity.node_intervals = node_intervals
    solverS.velocity.npts = VG.shape
    solverS.velocity.values = VG/VpVs
    return solverP.velocity, solverS.velocity


# @decoratortimer(2)
def TravelTimeTable(solver, src_idx, saveingPath):
    """Given velocity grid, it computes travel times by specifying source index

    Args:
        solver (pykonal.solver): a pykonal solver object
        src_idx (tuple): a tuple contains source location (indexed value not absolute)
        saveingPath (str): file for which the travel time table will be saved in
    """
    solver.traveltime.values[src_idx] = 0
    solver.unknown[src_idx] = False
    solver.trial.push(*src_idx)
    solver.solve()
    solver.traveltime.to_hdf(saveingPath)


# @decoratortimer(2)
def computeTTT(velocity, source_z, z_interval, saveingPath):
    """Creating travel times table

    Args:
        velocity (array): an array contains velocity grid
        source_z (float): source depth in km
        z_interval (float): node interval in Z direction
        saveingPath (str): file for which the travel time table will be saved in
    """
    zIndx = int(source_z/z_interval)
    src_idx = (0, 0, zIndx)
    solver = pykonal.EikonalSolver(coord_sys="cartesian")
    solver.velocity.min_coords = velocity.min_coords
    solver.velocity.node_intervals = velocity.node_intervals
    solver.velocity.npts = velocity.npts
    solver.velocity.values = velocity.values
    TravelTimeTable(solver, src_idx, saveingPath)

# @decoratortimer(2)


def generateTTT(velocityModelDict, velocityType, xGridMax, zGridMax, node_intervals=(1., 1., 1.), decimationFactor=0):
    """generate travel time tables

    Args:
        velocityModelDict (dict): a dictionary contains P velocity and Vp/Vs ration
        velocityType (str): type of velocity (P or S)
        xGridMax (int): maximum horizontal distance in km
        zGridMax (int): maximum vertical distance in km
        node_intervals (tuple, optional): node intervals in X,Y and Z direction. Defaults to (1., 1., 1.).
        decimationFactor (int, optional): decimation factor. Defaults to 0.
    """
    nXnYnZ = xGridMax, 1, zGridMax
    min_coords = (0.0, 0.0, 0.0)
    velocity_P, velocity_S = generateVelocityGrid(
        velocityModelDict, nXnYnZ, node_intervals, min_coords, decimationFactor)
    z_interval = node_intervals[-1]
    Path("ttt").mkdir(parents=True, exist_ok=True)
    for source_z in tzip(range(zGridMax)):
        source_z = source_z[0]
        saveingPath = os.path.join("ttt", "dep{z:003d}{vt:s}.hdf5".format(
            z=int(source_z), vt=velocityType))
        if os.path.exists(saveingPath):
            os.remove(saveingPath)
        if velocityType == "P":
            computeTTT(velocity_P, source_z, z_interval, saveingPath)
        else:
            computeTTT(velocity_S, source_z, z_interval, saveingPath)

# @decoratortimer(2)


def extractTT(traveltime, receivers=[]):
    """Extract travel time from given computational grid

    Args:
        traveltime (pykonal.solver.traveltime): on object contains travel time grid
        receivers (list, optional): a list contains station location. Defaults to [].

    Returns:
        _type_: _description_
    """
    extractedTT = list(map(traveltime.value, receivers))
    return extractedTT
