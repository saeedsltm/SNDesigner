from numpy import array, ones, max, float16
import pykonal
from math import ceil

def makeVelocityGrid(velocityModelDict, nXnYnZ, node_intervals=(1., 1., 1.), min_coords=(0. ,0. ,0.), df=0):
    """
    Make velocity grid for P and S velocities.
    - Inputs:
    vm: a dictionary containing 1D velocity model,
    nXnYnZ: number of grid points in x,y and Z directions,
    node_intervals: grid interval in x,y,z directions,
    min_coords: origin of the computational grid,
    df: decimation factor which is power of 2.
    - Outputs:
    solverP: P velocity grid,
    solverS: S velocity grid.
    """
    DATA_TYPE = float16
    vgP = ones((nXnYnZ), dtype=DATA_TYPE)
    maxZ = nXnYnZ[2] * node_intervals[2]
    Vp, Z, VpVs = [v for _,v in velocityModelDict.items()]
    for vp, z in zip(Vp, Z):
        indx = int(z * 1/node_intervals[2])
        if z >= maxZ:
            errMsg = "Velocity interface ({0}km) exceeds number of grids ({1}km) in Z direction!".format(z, maxZ)
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

def TravelTimeTable(solver, src_idx):
    """
    Given velocity grid, it computes travel times by specifying source index
    - Inputs:
    VG: 3D velocity grid,
    src_idx: source location based on computational grid index,
    min_coords: origin of the computational grid,
    node_intervals: intervals between computational grid points in x,y and z,
    df: decimation factor, which will be power of 2.
    - Outputs:
    travel time tables will be written to disk.
    """
    solver.traveltime.values[src_idx] = 0
    solver.unknown[src_idx] = False
    solver.trial.push(*src_idx)
    solver.solve()
    return solver.traveltime

def computeTTT(velocity, source_z, z_interval):
    """
    Creating travel times table.
    - Inputs:
    velocity: velocity object returned from pykonal,
    source_z: source depth in km,
    z_interval: grid interval in z direction,
    - Outputs:
    travelTime: pykonal traveltime will be returned
    """
    zIndx = int(source_z/z_interval)
    src_idx = (0 ,0 , zIndx)
    solver = pykonal.EikonalSolver(coord_sys="cartesian")
    solver.velocity.min_coords = velocity.min_coords
    solver.velocity.node_intervals = velocity.node_intervals
    solver.velocity.npts = velocity.npts
    solver.velocity.values = velocity.values
    tt = TravelTimeTable(solver, src_idx)
    return tt

def extractTT(traveltime, receivers=[]):
    """
    Extract travel time from given computational grid
    - Inputs:
    traveltime: travelTime object returned from pykonal
    receivers: list of desired station points (x,y,z) for getting travel times.
    - Outputs:
    tt: list of calculated travel times for desired station points.
    """
    extractedTT = list(map(traveltime.value, receivers))
    return extractedTT

def generateTT(velocityModelDict, velocityType, source_z, zGridMax, receivers, node_intervals, decimationFactor):
    nXnYnZ = ceil(max(receivers)+5), 1, zGridMax
    min_coords = (0.0, 0.0, 0.0)
    node_intervals = node_intervals
    decimationFactor = decimationFactor
    velocity_P, velocity_S = makeVelocityGrid(velocityModelDict, nXnYnZ, node_intervals, min_coords, decimationFactor)
    z_interval = node_intervals[-1]
    if velocityType == "P":
        traveltimeP = computeTTT(velocity_P, source_z, z_interval)
        extractedTT = extractTT(traveltimeP, receivers)
    else:
        traveltimeS = computeTTT(velocity_S, source_z, z_interval)
        extractedTT = extractTT(traveltimeS, receivers)
    return extractedTT