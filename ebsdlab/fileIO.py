# Read and write EBSD files: loaders fill an EBSD instance, writeANG writes one

from __future__ import annotations
import math
import os
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any
import numpy as np
from ._rotation import asBungeEulers
from .symmetry import Symmetry
if TYPE_CHECKING:
    from .ebsd import EBSD

# The low Laue classes m-3, 6/m, 4/m, -3 use the high ones m-3m, 6/mmm, 4/mmm, -3m of their crystal system,
#    which merges some distinct orientations; add separate lattices if such phases matter
TSL_SYMMETRIES = {43: 'cubic', 23: 'cubic', 62: 'hexagonal', 6: 'hexagonal', 42: 'tetragonal', 4: 'tetragonal',
                  22: 'orthorhombic', 32: 'trigonal', 3: 'trigonal', 2: 'monoclinic', 1: 'triclinic', 'm-3m': 'cubic'}
OXFORD_LAUE_GROUPS = {11: 'cubic', 10: 'cubic', 9: 'hexagonal', 8: 'hexagonal', 5: 'tetragonal', 4: 'tetragonal',
                      3: 'orthorhombic', 7: 'trigonal', 6: 'trigonal', 2: 'monoclinic', 1: 'triclinic'}


def loadANG(ebsd: EBSD, fileName: str = '') -> None:
    """Load .ang file: filename saved in ebsd. No need to use it

    Args:
       fileName: file to read [default: ebsd.fileName]
    """
    if fileName:
        ebsd.fileName = fileName
    print('Load .ang file: ', ebsd.fileName)
    keys = ['MaterialName', 'LatticeConstants', 'WorkingDistance', 'SEMVoltage', 'GRID:', 'Symmetry']
    fileHandle = open(ebsd.fileName)
    keyValues: list[Any] = [''] * len(keys)  # actual values
    symmetries: dict[int, Any] = {}  # phase number: TSL symmetry code
    phase = 0
    for line in fileHandle:
        if line[0:10] == '# OPERATOR':
            break
        parts = line.split()
        if parts[:2] == ['#', 'Phase'] and len(parts) == 3:  # phases may be listed in any order
            phase = int(parts[2])
        for key in keys:
            searchTerm = '# '+key
            if searchTerm == line[0:len(searchTerm)]:
                index = keys.index(key)
                value: Any = line.rstrip().split()[2:]
                if len(value) == 1:
                    value = value[0]
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                keyValues[index] = value
                if key == 'Symmetry':
                    symmetries[phase or len(symmetries)+1] = value
                break
    ebsd.meta = dict(list(zip(keys, keyValues)))
    ebsd.sym = [Symmetry(TSL_SYMMETRIES.get(symmetries.get(i, ''), ''))
                for i in range(1, max(symmetries, default=0)+1)]
    # read data
    data           = np.loadtxt(fileHandle)
    ebsd.phi1      = data[:, 0].astype(float)
    ebsd.phi       = data[:, 1].astype(float)
    ebsd.phi2      = data[:, 2].astype(float)
    ebsd.iq        = data[:, 5].astype(float)
    ebsd.ci        = data[:, 6].astype(float)
    ebsd.phaseID   = data[:, 7].astype(np.uint8)
    ebsd.phaseID  += not ebsd.phaseID.any()  # single-phase EDAX files use 0
    ebsd.semSignal = data[:, 8].astype(np.uint8)
    ebsd.fit       = data[:, 9].astype(float)
    # do coordinates
    ebsd._setGrid(data[:, 3], data[:, 4])
    ebsd.width     = max(data[:, 3])
    ebsd.height    = max(data[:, 4])
    ebsd.ratio     = ebsd.width/ebsd.height
    fileHandle.close()
    del data
    return


def loadTXT(ebsd: EBSD, fileName: str = '', update: bool = False) -> None:
    """ read txt file and possibly update data

    Args:
      fileName: fileName to load (partition data from OIM)
      update: update data or read new (read-new: default)
    """
    # TODO: Symmetry has to be read and used
    startTime = time.time()
    print('Load .txt file:', fileName)
    if not fileName:
        fileName = ebsd.fileName
    fileHandle = open(fileName)
    foundKeys: dict[str, int] = {}
    for line in fileHandle:
        if line[0] != '#':
            break
        parts = line.split()
        if len(parts) < 2:
            continue
        if parts[1] == 'Header:':
            print('   Header: ', parts[2])
            continue
        if 'Column' in parts:
            foundKeys[parts[3]] = int(parts[2].split(':')[0].split('-')[0])
    print('   Found data:', foundKeys)
    if 'Grain' in foundKeys:  # open new array if data exists
        ebsd.grainID = -np.ones_like(ebsd.phaseID)

    # read data
    data = np.loadtxt(fileName)
    print('   Reading file of size ', data.shape, '  this can take a bit...')
    if update:
        ebsd.mask[:] = False
        print('Warning: this is too slow')
        #  TODO: be intelligent where you seearch, check if old and new data monotonically increases then
        # search in sections of equal y or subdivide into half, of half of half
        xs, ys = ebsd.xy()
        for i in range(data.shape[0]):
            x = data[i, foundKeys['x,'] - 1].astype(float)
            y = data[i, foundKeys['x,'] - 0].astype(float)
            # identify index: nice and much much slowes
            idx = np.argmax(np.logical_and(np.abs(xs-x) < ebsd.stepSizeX/10.0, # very safe error of STEPSIZE/10
                                           np.abs(ys-y) < ebsd.stepSizeX/10.0))
            # update
            ebsd.mask[idx] = True
            ebsd.phi1[idx] = data[i, foundKeys['phi1,'] - 1].astype(np.float16)
            ebsd.phi[idx]  = data[i, foundKeys['phi1,'] - 0].astype(np.float16)
            ebsd.phi2[idx] = data[i, foundKeys['phi1,'] + 1].astype(np.float16)
            if 'IQ' in foundKeys:
                ebsd.iq[idx] = data[i, foundKeys['IQ']  - 1].astype(np.float16)
            if 'CI' in foundKeys:
                ebsd.ci[idx] = data[i, foundKeys['CI']  - 1].astype(np.float16)
            if 'Fit' in foundKeys:
                ebsd.fit[idx] = data[i, foundKeys['Fit'] - 1].astype(np.float16)
            if 'Phase' in foundKeys:
                ebsd.phaseID[idx] = data[i, foundKeys['Phase'] - 1].astype(np.float16)
            if 'sem' in foundKeys:
                ebsd.semSignal[idx] = data[i, foundKeys['sem'] - 1].astype(np.float16)
            if 'Grain' in foundKeys:
                ebsd.grainID[idx] = data[i, foundKeys['Grain'] - 1].astype(np.float16)
            # stepSizeX, width, height etc do not change
    else:  # read new
        ebsd.phi1 = data[:, foundKeys['phi1,'] - 1].astype(np.float16)
        ebsd.phi = data[:, foundKeys['phi1,'] - 0].astype(np.float16)
        ebsd.phi2 = data[:, foundKeys['phi1,'] + 1].astype(np.float16)
        x = data[:, foundKeys['x,'] - 1].astype(float)
        y = data[:, foundKeys['x,'] - 0].astype(float)
        if 'IQ' in foundKeys:
            ebsd.iq = data[:, foundKeys['IQ'] - 1].astype(np.float16)
        if 'CI' in foundKeys:
            ebsd.ci = data[:, foundKeys['CI'] - 1].astype(np.float16)
        if 'Fit' in foundKeys:
            ebsd.fit = data[:, foundKeys['Fit'] - 1].astype(np.float16)
        if 'Phase' in foundKeys:
            ebsd.phaseID = data[:,foundKeys['Phase'] - 1].astype(np.uint8)
            ebsd.phaseID += not ebsd.phaseID.any()  # single-phase EDAX files use 0
        if 'sem' in foundKeys:
            ebsd.semSignal = data[:, foundKeys['sem'] - 1].astype(np.float16)
        ebsd.mask   = np.ones_like(x, dtype=bool)
        ebsd.width  = max(x)
        ebsd.height = max(y)
        ebsd.ratio  = ebsd.width/ebsd.height
        ebsd._setGrid(x, y)
    fileHandle.close()
    print('Duration loadTXT: ', int(np.round(time.time()-startTime)), 'sec')
    return


def loadOSC(ebsd: EBSD, fileName: str = '') -> None:
    """Load .osc file. Copied from mtex and translated into python
    TODO: SEMsignal not parsed correctly

    Args:
       fileName: file to read [default: ebsd.fileName]
    """
    print('TODO: Symmetry has to be read and used')
    if fileName:
        ebsd.fileName = fileName
    print('Load .osc file: ', ebsd.fileName)

    # OSC stores its numeric values as little-endian 32-bit values.
    # Using NumPy's platform-sized ``float`` (normally float64) desynchronizes the reader after the first
    #    step-size field.
    startBytes = bytes.fromhex('B9 0B EF FF 02 00 00 00')
    raw        = Path(ebsd.fileName).read_bytes()
    header     = np.frombuffer(raw, dtype='<u4', count=8)
    n          = int(header[6])  # number of data points
    startPosition = raw.find(startBytes)
    if startPosition < 0:
        raise ValueError('OSC data-block marker was not found.')

    # OSC versions differ: a uint32 count field may precede the x and y step sizes, records have 10 or more
    # columns (e.g. PRIAS or EDS), and other blocks may follow. As in mtex, the layout is found where the
    # second record lies one step along x: x = stepSizeX, y = 0.
    base = startPosition + len(startBytes)
    for dataOffset, nColumns in ((offset, columns) for offset in (base, base+4) for columns in range(10, 31)):
        if dataOffset + 8 + n*nColumns*4 > len(raw):
            continue
        stepX, stepY = np.frombuffer(raw, dtype='<f4', count=2, offset=dataOffset).astype(float)
        second = np.frombuffer(raw, dtype='<f4', count=5, offset=dataOffset + 8 + nColumns*4)
        if stepX > 1e-6 and np.isclose(second[3], stepX, rtol=1e-4) and second[4] == 0:
            break
    else:
        raise ValueError('OSC data layout was not recognized.')
    ebsd.stepSizeX, ebsd.stepSizeY = stepX, stepY
    data = np.frombuffer(raw, dtype='<f4', count=n*nColumns, offset=dataOffset + 8).reshape(n, nColumns)
    ebsd.phi1 = data[:, 0].astype(np.float16)
    ebsd.phi  = data[:, 1].astype(np.float16)
    ebsd.phi2 = data[:, 2].astype(np.float16)
    # scale image-quality to float16
    iqScale   = 10.0**max(0, math.ceil(math.log10(max(float(data[:, 5].max()), 1.0)/float(np.finfo(np.float16).max))))
    if iqScale > 1:
        ebsd.iq = (data[:, 5]/iqScale).astype(np.float16)
    ebsd.ci        = data[:, 6].astype(np.float16)
    ebsd.phaseID   = data[:, 7].astype(np.uint8)
    ebsd.phaseID  += not ebsd.phaseID.any()  # single-phase EDAX files use 0
    ebsd.semSignal = data[:, 8].astype(np.float16)  # SEMSignal
    ebsd.fit       = data[:, 9].astype(np.float16)  # Fit
    ebsd.width     = float(max(data[:, 3]))
    ebsd.height    = float(max(data[:, 4]))
    ebsd.ratio     = ebsd.width/ebsd.height
    ebsd._setGrid(data[:, 3].astype(float), data[:, 4].astype(float))
    del data
    return


def loadCRC(ebsd: EBSD, fileName: str = '') -> None:
    """Load .crc file; filename saved in ebsd. No need to use it. Copied from mtex and translated into python

    Args:
       fileName: file to read [default: ebsd.fileName]; the .cpr file of the same name holds the metadata
    """
    if fileName:
        ebsd.fileName = fileName
    cprFileName = ebsd.fileName[:-4]+'.cpr'
    print('Load .crc file: ', ebsd.fileName, cprFileName)
    if not os.path.exists(cprFileName):
        print('CPR file does not exist')
    cprFile = open(cprFileName)
    cprData: dict[str, dict[str, Any]] = {}
    for line in cprFile:
        line = line.strip()
        if not line:
            continue
        if line[0] == '[':
            title = line[1:-1].lower()
            cprData[title] = {}
            continue
        key, value = line.split('=')[0], line.split('=')[1]
        try:
            cprData[title][key.lower()] = float(value)
        except ValueError:
            cprData[title][key.lower()] = value.lower()
    cprFile.close()
    # print "META DATA",cprData
    if 'griddistx' not in cprData['job']:
        raise ValueError(f'CRC file is not a grid map (JobMode={cprData["general"].get("jobmode")})')
    ebsd.stepSizeX = np.double(cprData['job']['griddistx'])
    ebsd.stepSizeY = np.double(cprData['job']['griddisty'])
    xcells         = int(cprData['job']['xcells'])
    ycells         = int(cprData['job']['ycells'])
    numDataPoints  = xcells * ycells
    ebsd.width     = xcells * ebsd.stepSizeX
    ebsd.height    = ycells * ebsd.stepSizeY
    ebsd.ratio     = ebsd.width/ebsd.height
    ebsd.sym = [Symmetry(OXFORD_LAUE_GROUPS.get(cprData[f'phase{i}']['lauegroup'], ''))
                for i in range(1, int(cprData['phases']['count'])+1)]
    # verify that data in correct order
    allColumnNames = [
        'X',                  # 1    4 bytes
        'Y',                  # 2       "
        'phi1',               # 3       "
        'Phi',                # 4       "
        'phi2',               # 5       "
        'MAD',                # 6       "
        'BC',                 # 7    1 byte
        'BS',                 # 8       "
        'Unknown',            # 9       "
        'Bands',              # 10      "
        'Error',              # 11      "
        'ReliabilityIndex']    # 12      "
    allDataType     = np.ones((12,), dtype=int)
    allDataType[:6] = 4
    allDataType[-1] = 4
    columnNames, columnType = ['Phase'], [1]
    for k in range(int(cprData['fields']['count'])):
        order = int(cprData['fields']['field'+str(k+1)])-1
        if order <= 12:
            columnNames.append(allColumnNames[order])
            columnType.append(allDataType[order])
        else:
            columnNames.append('Unknown'+str(order))
            columnType.append(4)
    expectedColumns = [ 'Phase', 'phi1', 'Phi', 'phi2', 'MAD', 'BC', 'BS', 'Bands', 'Error', 'ReliabilityIndex']
    if columnNames == expectedColumns:
        print('  CRC-Data in correct order')
    else:
        print('  WARNING! CRC-Data not in correct order! WARNING')
        print(f'    should be {expectedColumns}')
        print('    is       ', columnNames)
        print('    if data missing at end, no problem')
    # coordinates
    xCoordinates = np.arange(xcells)*ebsd.stepSizeX
    yCoordinates = np.arange(ycells)*ebsd.stepSizeY
    x, y         = np.meshgrid(xCoordinates, yCoordinates)
    ebsd._setGrid(x.flatten(), y.flatten())

    # read data from crcFile: packed little-endian records
    recordType = [('phase', 'u1'), ('phi1', '<f4'), ('phi', '<f4'), ('phi2', '<f4'), ('ci', '<f4'),
                  ('bc', 'u1'), ('bs', 'u1'), ('bands', 'u1'), ('error', 'u1')]
    if 'ReliabilityIndex' in columnNames:
        recordType.append(('ri', '<f4'))
    recordType += [(name, '<f4') for name in columnNames if name.startswith('Unknown') and name != 'Unknown']
    data         = np.fromfile(ebsd.fileName, dtype=recordType, count=numDataPoints)
    ebsd.phaseID = data['phase'].copy()
    ebsd.bc, ebsd.bs, ebsd.bands, ebsd.error = (data[i].copy() for i in ('bc', 'bs', 'bands', 'error'))
    ebsd.phi1, ebsd.phi, ebsd.phi2, ebsd.ci  = (data[i].astype(float) for i in ('phi1', 'phi', 'phi2', 'ci'))
    ebsd.ri = data['ri'].astype(float) if 'ri' in (data.dtype.names or ()) else np.zeros(numDataPoints)
    ebsd.iq, ebsd.semSignal, ebsd.fit        = (np.zeros(numDataPoints) for _ in range(3))
    if np.max(ebsd.phaseID) > len(ebsd.sym):
        print('ERROR in reading CRC: symmetries do not match', len(ebsd.sym), np.max(ebsd.phaseID))
    return


def loadCTF(ebsd: EBSD, fileName: str = '') -> None:
    """Load Oxford .ctf text file; filename saved in ebsd

    Args:
       fileName: file to read [default: ebsd.fileName]
    """
    if fileName:
        ebsd.fileName = fileName
    print('Load .ctf file: ', ebsd.fileName)
    # some AZtec versions write decimal commas; ',' is no field separator in .ctf
    with open(ebsd.fileName, encoding='utf-8-sig', errors='replace') as fileHandle:
        lines = fileHandle.read().replace(',', '.').splitlines()
    header: dict[str, str] = {}
    for iLine, line in enumerate(lines):
        parts = line.split('\t')
        if parts[0] == 'Phases':
            ebsd.sym = [Symmetry(OXFORD_LAUE_GROUPS.get(int(phaseLine.split('\t')[3]), ''))
                        for phaseLine in lines[iLine+1:iLine+1+int(parts[1])]]
        elif parts[0] == 'Phase':  # column names, data follows
            break
        elif len(parts) > 1:
            header[parts[0]] = parts[1]
    if header.get('JobMode') != 'Grid':
        raise ValueError(f'CTF file is not a grid map (JobMode={header.get("JobMode")}); only grid maps allowed')
    ebsd.stepSizeX, ebsd.stepSizeY = float(header['XStep']), float(header['YStep'])
    xcells, ycells                 = int(header['XCells']), int(header['YCells'])
    data                           = np.loadtxt(lines[iLine+1:], ndmin=2)
    if len(data) != xcells*ycells:
        raise ValueError(f'CTF file has {len(data)} points, but XCells*YCells = {xcells*ycells}.')
    ebsd.width, ebsd.height        = xcells*ebsd.stepSizeX, ycells*ebsd.stepSizeY
    ebsd.ratio                     = ebsd.width/ebsd.height
    # coordinates from the header: the X and Y columns are rounded
    x, y = np.meshgrid(np.arange(xcells)*ebsd.stepSizeX, np.arange(ycells)*ebsd.stepSizeY)
    ebsd._setGrid(x.flatten(), y.flatten())
    # columns: Phase X Y Bands Error Euler1 Euler2 Euler3 MAD BC BS
    ebsd.phaseID                   = data[:, 0].astype(np.uint8)
    ebsd.bands, ebsd.error         = data[:, 3].astype(np.uint8), data[:, 4].astype(np.uint8)
    ebsd.phi1, ebsd.phi, ebsd.phi2 = np.radians(data[:, 5:8]).T
    ebsd.ci                        = data[:, 8]
    ebsd.bc, ebsd.bs               = data[:, 9].astype(np.uint8), data[:, 10].astype(np.uint8)
    ebsd.ri, ebsd.iq, ebsd.semSignal, ebsd.fit = (np.zeros(len(data)) for _ in range(4))
    return


def loadVoid(ebsd: EBSD, rotation: str) -> None:
    """rotation angles in degree

    Args:
       rotation: Euler angles "phi1|Phi|phi2" in degrees, optionally followed by "|spread" (standard deviation
                 of the random scatter in radians) and "|numberPerAxis" [default: 6]; without "|": no rotation
    """
    numPerAxis, distrib = 6, 0.0
    if '|' in rotation:
        values = [float(i) for i in rotation.split('|')]
        if len(values) == 3:
            phi1, phi, phi2 = np.radians(values)
        elif len(values) == 4:
            phi1, phi, phi2 = np.radians(values[:3])
            distrib = values[-1]
        elif len(values) == 5:
            phi1, phi, phi2 = np.radians(values[:3])
            distrib, numPerAxis = values[-2], int(values[-1])
        else:
            print('ERROR')
            return
        print('   Euler angles:', np.round(phi1, 2), np.round(phi, 2), np.round(phi2, 2),
              '| distribution:', distrib, '| numberPerAxis:', numPerAxis)
    else:
        phi1, phi, phi2 = 0, 0, 0
    if distrib < 0.001:
        distrib = 0.001
    ebsd.sym.append(Symmetry('cubic'))
    ebsd.stepSizeX = 1.
    numDataPoints  = int(numPerAxis**2)
    coordinates    = np.arange(numPerAxis)*ebsd.stepSizeX
    x, y           = np.meshgrid(coordinates, coordinates)
    ebsd._setGrid(x.flatten(), y.flatten())
    ebsd.phaseID   = np.ones((numDataPoints), dtype=np.uint8)
    ebsd.phi1      = np.zeros((numDataPoints), dtype=float)+phi1 + \
                        np.random.normal(loc=0, scale=distrib, size=numDataPoints)
    ebsd.phi       = np.zeros((numDataPoints), dtype=float)+phi + \
                        np.random.normal(loc=0, scale=distrib, size=numDataPoints)
    ebsd.phi2      = np.zeros((numDataPoints), dtype=float)+phi2 + \
                        np.random.normal(loc=0, scale=distrib, size=numDataPoints)
    ebsd.ci        = np.ones((numDataPoints), dtype=float)
    ebsd.width     = np.max(x)
    ebsd.height    = np.max(y)
    ebsd.ratio     = ebsd.width/ebsd.height
    return


def writeANG(ebsd: EBSD, fileName: str) -> None:
    """write body of ang file

    Args:
       fileName: file name
    """
    startTime = time.time()
    fileOut = open(fileName, 'w')
    fileOut.write('# MaterialName void\n')
    fileOut.write('# Formula \n')
    # adopt for HCP (fcc and bcc the same)
    fileOut.write('# Symmetry 43\n')
    fileOut.write('# LatticeConstants 1.0 1.0 1.0 90.0 90.0 90.0\n')
    fileOut.write('# NumberFamilies 4\n')
    fileOut.write('# khlFamilies 1 1 1 1 0.0\n')  # adopt for HCP
    fileOut.write('# khlFamilies 2 0 0 1 0.0\n')
    fileOut.write('# khlFamilies 2 2 0 1 0.0\n')
    fileOut.write('# khlFamilies 3 1 1 1 0.0\n')
    fileOut.write(f'#\n# GRID: {ebsd.grid}\n#\n')
    xs, ys = ebsd.xy()
    phaseID = ebsd.phaseID if ebsd.phaseID.max() > 1 else np.zeros_like(ebsd.phaseID)
    for i in range(ebsd.nPoints):
        phi1, phi, phi2 = tuple(asBungeEulers(ebsd.quaternions[i]))
        fileOut.write(f' {phi1:8.5f} {phi:8.5f} {phi2:8.5f} {xs[i]:12.5f} {ys[i]:12.5f} {ebsd.iq[i]:8.3f}'
                      f' {ebsd.ci[i]:6.3f} {phaseID[i]:2d} {ebsd.semSignal[i]:6d} {ebsd.fit[i]:7.3f}\n')
    fileOut.close()
    print('Duration writeANG: ', int(np.round(time.time()-startTime)), 'sec')
    return

LOADERS = {'.ang': loadANG, '.osc': loadOSC, '.txt': loadTXT, '.crc': loadCRC, '.ctf': loadCTF}
