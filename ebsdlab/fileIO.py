"""Read and write EBSD files: loaders fill an EBSD instance, writeANG writes one

Conventions (sample frame, crystal frame, vendor frames): docs/source/conventions.rst
"""

from __future__ import annotations
import itertools
import math
import os
import re
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any
import numpy as np
from scipy.spatial.transform import Rotation
from ._rotation import asBungeEulers, fromBungeEulers, multiply
from .symmetry import GROUPS, Symmetry
if TYPE_CHECKING:
    from .ebsd import EBSD

# vendor frames, see docs/source/conventions.rst: q = SAMPLE * q_file * CRYSTAL
EDAX_SAMPLE    = Rotation.from_rotvec(np.pi*np.array([1, -1, 0])/np.sqrt(2))  # 180° about [1-10]
OXFORD_CRYSTAL = Rotation.from_euler('z', -30, degrees=True)                  # hexagonal/trigonal only
FRAMES = {'.ang': (EDAX_SAMPLE, None), '.osc': (EDAX_SAMPLE, None), '.txt': (EDAX_SAMPLE, None),
          '.h5': (EDAX_SAMPLE, None), '.crc': (None, OXFORD_CRYSTAL), '.ctf': (None, OXFORD_CRYSTAL)}


def rotateToConventions(ebsd: EBSD, suffix: str) -> None:
    """Rotate the orientations, as stored in a file, into the conventions of docs/source/conventions.rst

    Args:
       ebsd: instance with quaternions and symmetries of the loaded file
       suffix: file suffix, e.g. '.ang'; files not in FRAMES are not rotated
    """
    sample, crystal = FRAMES.get(suffix, (None, None))
    if sample is not None:
        ebsd.quaternions = multiply(sample, ebsd.quaternions)
    if crystal is not None:
        for phase, sym in enumerate(ebsd.sym):
            points = ebsd.phaseID == phase
            if sym.lattice in ('hexagonal', 'trigonal') and points.any():
                ebsd.quaternions[points] = multiply(ebsd.quaternions[points], crystal)

# The low Laue classes m-3, 6/m, 4/m, -3 use the high ones m-3m, 6/mmm, 4/mmm, -3m of their crystal system;
#    the loaders warn when they meet one.
# - The result: the extra operations of the high class merge distinct orientations, so misorientations and KAM
#   of such phases are underestimated and IPF colors do not separate those orientations.
# - Who has low Laue classes: mostly minerals and oxides
#    - m-3 (point groups 23, m-3): pyrite FeS2, bixbyite (Mn2O3, Y2O3, In2O3/ITO), skutterudite CoAs3,
#      ullmannite, NaClO3, langbeinite
#    - 6/m (6, -6, 6/m): apatite incl. hydroxyapatite (bone, teeth), nepheline, cancrinite, UCl3-type
#    - 4/m (4, -4, 4/m): scheelite CaWO4, wulfenite, stolzite, LiYF4, scapolite, leucite;
#      intermetallics Ni4Mo, Ni4W (D1a)
#    - -3 (3, -3): dolomite, ilmenite FeTiO3, geikielite MgTiO3, willemite, phenakite, dioptase
# - Full support needs the rotation groups T, C6, C4, C3 and their standard triangles in symmetry.GROUPS.
TSL_SYMMETRIES = {43: 'cubic', 23: 'cubic', 62: 'hexagonal', 6: 'hexagonal', 42: 'tetragonal', 4: 'tetragonal',
                  22: 'orthorhombic', 32: 'trigonal', 3: 'trigonal', 2: 'monoclinic', 1: 'triclinic', 'm-3m': 'cubic'}
# writeANG: first code of each lattice, i.e. the high Laue class
TSL_CODES = {lattice: code for code, lattice in reversed(TSL_SYMMETRIES.items())}
OXFORD_LAUE_GROUPS = {11: 'cubic', 10: 'cubic', 9: 'hexagonal', 8: 'hexagonal', 5: 'tetragonal', 4: 'tetragonal',
                      3: 'orthorhombic', 7: 'trigonal', 6: 'trigonal', 2: 'monoclinic', 1: 'triclinic'}
LOW_LAUE_CLASSES = {'TSL': {23: 'm-3', 6: '6/m', 4: '4/m', 3: '-3'}, 'Oxford': {10: 'm-3', 8: '6/m', 4: '4/m', 6: '-3'}}


def symmetryFromCode(code: Any, convention: str) -> Symmetry:
    """Symmetry of a phase from its file code; warns for low Laue classes, which use the high one

    Args:
       code: symmetry code in the file
       convention: 'TSL' (.ang, .txt, .osc) or 'Oxford' (.crc, .ctf Laue group)

    Returns:
       symmetry; not identified if the code is unknown
    """
    if code in LOW_LAUE_CLASSES[convention]:
        print(f'   Warning: Laue class {LOW_LAUE_CLASSES[convention][code]} is treated as the higher class of its '
              'crystal system; see fileIO.TSL_SYMMETRIES')
    return Symmetry((TSL_SYMMETRIES if convention == 'TSL' else OXFORD_LAUE_GROUPS).get(code, ''))


def tslSymmetries(headerLines: list[str]) -> list[Symmetry]:
    """Phase symmetries from a TSL header (.ang, .txt): '# Phase k' is followed by '# Symmetry code'

    Args:
       headerLines: header lines of the file

    Returns:
       symmetry of phase 1, 2, ...; phases without symmetry are not identified
    """
    symmetries: dict[int, Any] = {}  # phase number: TSL symmetry code
    phase = 0
    for line in headerLines:
        parts = line.split()
        if parts[:2] == ['#', 'Phase'] and len(parts) == 3:  # phases may be listed in any order
            phase = int(parts[2])
        elif parts[:2] == ['#', 'Symmetry'] and len(parts) > 2:
            symmetries[phase or len(symmetries)+1] = int(parts[2]) if parts[2].isdigit() else parts[2]
    return [symmetryFromCode(symmetries.get(i, ''), 'TSL') for i in range(1, max(symmetries, default=0)+1)]


def oscSymmetries(header: bytes) -> list[Symmetry]:
    """Phase symmetries from the header block of an .osc file. As in mtex, the block is not self describing:
    a phase record is where a 256-byte name is followed by a TSL symmetry code (int32) and six plausible
    lattice constants (float32).

    Args:
       header: bytes between the header marker and the data-block marker

    Returns:
       symmetry of phase 1, 2, ... in the order of the file
    """
    symmetries, pos = [], 0
    while pos + 288 <= len(header):
        name = header[pos:pos+256].split(b'\0')[0]
        code = int(np.frombuffer(header, '<i4', 1, pos+256)[0])
        cell = np.frombuffer(header, '<f4', 6, pos+260)
        if re.fullmatch(rb'[ -~]*[A-Za-z][ -~]*', name) and 0 < code <= 131 and \
                np.all((cell[:3] > 0) & (cell[:3] < 1000)) and np.all((cell[3:] > 0) & (cell[3:] < 180)):
            symmetries.append(symmetryFromCode(code, 'TSL'))
            pos += 288
        else:
            pos += 1
    return symmetries


def loadANG(ebsd: EBSD, fileName: str = '') -> None:
    """Load .ang file: filename saved in ebsd. No need to use it

    Args:
       fileName: file to read [default: ebsd.fileName]
    """
    if fileName:
        ebsd.fileName = fileName
    print('Load .ang file: ', ebsd.fileName)
    keys = ['MaterialName', 'LatticeConstants', 'WorkingDistance', 'SEMVoltage', 'GRID:', 'Symmetry']
    with open(ebsd.fileName, encoding='utf-8', errors='replace') as fileHandle:  # header; loadtxt reads the data
        headerLines = list(itertools.takewhile(lambda line: line.startswith('#'), fileHandle))
    keyValues: list[Any] = [''] * len(keys)  # actual values
    for line in headerLines:
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
                break
    ebsd.meta = dict(list(zip(keys, keyValues)))
    ebsd.sym  = tslSymmetries(headerLines)
    # read data
    data           = np.loadtxt(ebsd.fileName)
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
    ebsd.setGrid(data[:, 3], data[:, 4])
    ebsd.width     = max(data[:, 3])
    ebsd.height    = max(data[:, 4])
    del data
    return


def loadH5(ebsd: EBSD, fileName: str = '') -> None:
    """Load EDAX OIM .h5 file, the binary counterpart of .ang; the first scan of the file is read.
    Needs h5py: pip install 'ebsdlab[h5]'

    Args:
       fileName: file to read [default: ebsd.fileName]
    """
    try:
        import h5py  # pylint: disable=import-outside-toplevel
    except ImportError as error:
        raise ImportError("Reading .h5 files needs h5py: pip install 'ebsdlab[h5]'") from error
    if fileName:
        ebsd.fileName = fileName
    print('Load .h5 file: ', ebsd.fileName)
    with h5py.File(ebsd.fileName, 'r') as fileHandle:
        manufacturer = fileHandle[' Manufacturer'][0].decode() if ' Manufacturer' in fileHandle else 'unknown'
        if manufacturer != 'EDAX':
            raise ValueError(f'Only EDAX OIM .h5 files are supported, not {manufacturer}')
        scans = [name for name in fileHandle if isinstance(fileHandle[name], h5py.Group) and 'EBSD' in fileHandle[name]]
        if len(scans) > 1:
            print('   File has several scans, read the first:', scans)
        header, data = fileHandle[scans[0]]['EBSD/Header'], fileHandle[scans[0]]['EBSD/Data']
        ebsd.meta = {'Scan': scans[0], 'Version': fileHandle[' Version'][0].decode(),
                     'GRID:': header['Grid Type'][0].decode()}
        phases = header['Phase']
        ebsd.sym = [symmetryFromCode(int(phases[k]['Symmetry'][0]), 'TSL')
                    for k in sorted(phases, key=int)]
        x, y = data['X Position'][()], data['Y Position'][()]
        valid = (x > -1e6) & (y > -1e6)  # invalid points: x = y = -1111111
        ebsd.phi1      = data['Phi1'][()][valid].astype(float)
        ebsd.phi       = data['Phi'][()][valid].astype(float)
        ebsd.phi2      = data['Phi2'][()][valid].astype(float)
        ebsd.iq        = data['IQ'][()][valid].astype(float)
        ebsd.ci        = data['CI'][()][valid].astype(float)
        ebsd.phaseID   = data['Phase'][()][valid].astype(np.uint8)
        ebsd.phaseID  += not ebsd.phaseID.any()  # single-phase EDAX files use 0
        ebsd.semSignal = data['SEM Signal'][()][valid].astype(np.uint8)
        ebsd.fit       = data['Fit'][()][valid].astype(float)
    if not valid.all():
        print('   Invalid points skipped:', (~valid).sum())
    ebsd.setGrid(x[valid].astype(float), y[valid].astype(float))
    ebsd.width     = float(x[valid].max())
    ebsd.height    = float(y[valid].max())
    return


def loadTXT(ebsd: EBSD, fileName: str = '', update: bool = False) -> None:
    """ read txt file and possibly update data

    Args:
      fileName: fileName to load (partition data from OIM)
      update: update data or read new (read-new: default)
    """
    startTime = time.time()
    print('Load .txt file:', fileName)
    if not fileName:
        fileName = ebsd.fileName
    with open(fileName, encoding='utf-8', errors='replace') as fileHandle:
        headerLines = list(itertools.takewhile(lambda line: line.startswith('#'), fileHandle))
    foundKeys: dict[str, int] = {}
    for line in headerLines:
        parts = line.split()
        if len(parts) < 2:
            continue
        if parts[1] == 'Header:':
            print('   Header: ', parts[2])
            continue
        if 'Column' in parts:
            foundKeys[parts[3]] = int(parts[2].split(':')[0].split('-')[0])
    print('   Found data:', foundKeys)
    if not update:
        ebsd.sym = tslSymmetries(headerLines)

    # read data
    data = np.loadtxt(fileName, ndmin=2)
    print('   Reading file of size ', data.shape, '  this can take a bit...')
    if update:  # e.g. a partition of the loaded map: only its points stay visible
        x, y   = data[:, foundKeys['x,'] - 1], data[:, foundKeys['x,'] - 0]
        idx    = ebsd.nearestIndex(x, y)
        xs, ys = ebsd.xy(idx)
        onMap  = (np.abs(xs-x) < ebsd.stepSizeX/10) & (np.abs(ys-y) < ebsd.stepSizeY/10)
        if not onMap.all():
            print('   Warning: points not on the loaded map are skipped:', (~onMap).sum())
        idx, data = idx[onMap], data[onMap]
        ebsd.mask[:]   = False
        ebsd.mask[idx] = True
        eulers = data[:, foundKeys['phi1,']-1:foundKeys['phi1,']+2]
        ebsd.quaternions[idx] = multiply(EDAX_SAMPLE, fromBungeEulers(eulers))
        if 'IQ' in foundKeys:
            ebsd.iq[idx] = data[:, foundKeys['IQ'] - 1]
        if 'CI' in foundKeys:
            ebsd.ci[idx] = data[:, foundKeys['CI'] - 1]
        if 'Fit' in foundKeys:
            ebsd.fit[idx] = data[:, foundKeys['Fit'] - 1]
        if 'Phase' in foundKeys:
            phaseID = data[:, foundKeys['Phase'] - 1].astype(np.uint8)
            ebsd.phaseID[idx] = phaseID + (not phaseID.any())  # single-phase EDAX files use 0
        if 'sem' in foundKeys:
            ebsd.semSignal[idx] = data[:, foundKeys['sem'] - 1]
        if 'Grain' in foundKeys:
            if not ebsd.grainID.size:
                ebsd.grainID = np.full(ebsd.nPoints, -1)
            ebsd.grainID[idx] = data[:, foundKeys['Grain'] - 1]
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
        if 'Grain' in foundKeys:
            ebsd.grainID = data[:, foundKeys['Grain'] - 1].astype(int)
        ebsd.mask   = np.ones_like(x, dtype=bool)
        ebsd.width  = max(x)
        ebsd.height = max(y)
        ebsd.setGrid(x, y)
    fileHandle.close()
    print('Duration loadTXT: ', int(np.round(time.time()-startTime)), 'sec')
    return


def loadOSC(ebsd: EBSD, fileName: str = '') -> None:
    """Load .osc file. Copied from mtex and translated into python
    TODO: SEMsignal not parsed correctly

    Args:
       fileName: file to read [default: ebsd.fileName]
    """
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
    headerStart = raw.find(bytes.fromhex('B9 0B EF FF 01 00 00 00'))
    ebsd.sym    = oscSymmetries(raw[headerStart+8:startPosition]) if 0 <= headerStart < startPosition else []

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
    # scale image-quality to float16; by design, the scaled IQ is used everywhere, also by writeANG
    float16Max = float(np.finfo(np.float16).max)  # pylint: disable=no-member
    iqScale   = 10.0**max(0, math.ceil(math.log10(max(float(data[:, 5].max()), 1.0)/float16Max)))
    ebsd.iq   = (data[:, 5]/iqScale).astype(np.float16)
    ebsd.meta['iqScale'] = iqScale
    ebsd.ci        = data[:, 6].astype(np.float16)
    ebsd.phaseID   = data[:, 7].astype(np.uint8)
    ebsd.phaseID  += not ebsd.phaseID.any()  # single-phase EDAX files use 0
    ebsd.semSignal = data[:, 8].astype(np.float16)  # SEMSignal
    ebsd.fit       = data[:, 9].astype(np.float16)  # Fit
    ebsd.width     = float(max(data[:, 3]))
    ebsd.height    = float(max(data[:, 4]))
    ebsd.setGrid(data[:, 3].astype(float), data[:, 4].astype(float))
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
    with open(cprFileName, encoding='utf-8', errors='replace') as cprFile:
        cprLines = cprFile.read().splitlines()
    cprData: dict[str, dict[str, Any]] = {}
    for line in cprLines:
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
    # print "META DATA",cprData
    if 'griddistx' not in cprData['job']:
        raise ValueError(f'CRC file is not a grid map (JobMode={cprData["general"].get("jobmode")}); '
                         'only grid maps are supported')
    ebsd.stepSizeX = np.double(cprData['job']['griddistx'])
    ebsd.stepSizeY = np.double(cprData['job']['griddisty'])
    xcells         = int(cprData['job']['xcells'])
    ycells         = int(cprData['job']['ycells'])
    numDataPoints  = xcells * ycells
    ebsd.width     = xcells * ebsd.stepSizeX
    ebsd.height    = ycells * ebsd.stepSizeY
    ebsd.sym = [symmetryFromCode(cprData[f'phase{i}']['lauegroup'], 'Oxford')
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
    ebsd.setGrid(x.flatten(), y.flatten())

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
            ebsd.sym = [symmetryFromCode(int(phaseLine.split('\t')[3]), 'Oxford')
                        for phaseLine in lines[iLine+1:iLine+1+int(parts[1])]]
        elif parts[0] == 'Phase':  # column names, data follows
            break
        elif len(parts) > 1:
            header[parts[0]] = parts[1]
    else:
        raise ValueError('CTF file has no data: the column line "Phase X Y ..." is missing')
    if header.get('JobMode') != 'Grid':
        raise ValueError(f'CTF file is not a grid map (JobMode={header.get("JobMode")}); only grid maps are supported')
    ebsd.stepSizeX, ebsd.stepSizeY = float(header['XStep']), float(header['YStep'])
    xcells, ycells                 = int(header['XCells']), int(header['YCells'])
    data                           = np.loadtxt(lines[iLine+1:], ndmin=2)
    if len(data) != xcells*ycells:
        raise ValueError(f'CTF file has {len(data)} points, but XCells*YCells = {xcells*ycells}.')
    ebsd.width, ebsd.height        = xcells*ebsd.stepSizeX, ycells*ebsd.stepSizeY
    # coordinates from the header: the X and Y columns are rounded
    x, y = np.meshgrid(np.arange(xcells)*ebsd.stepSizeX, np.arange(ycells)*ebsd.stepSizeY)
    ebsd.setGrid(x.flatten(), y.flatten())
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
    distrib = max(distrib, 0.001)
    ebsd.sym.append(Symmetry('cubic'))
    ebsd.stepSizeX = 1.
    numDataPoints  = int(numPerAxis**2)
    coordinates    = np.arange(numPerAxis)*ebsd.stepSizeX
    x, y           = np.meshgrid(coordinates, coordinates)
    ebsd.setGrid(x.flatten(), y.flatten())
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
    return


def writeANG(ebsd: EBSD, fileName: str) -> None:
    """write body of ang file

    Args:
       fileName: file name
    """
    startTime = time.time()
    with open(fileName, 'w', encoding='utf-8') as fileOut:
        # one block per phase: TSL code of the crystal system (0: not identified) and its default unit cell
        for phase, sym in enumerate(ebsd.sym[1:], start=1):
            cell      = GROUPS[sym.lattice]['cell'] if sym.lattice else {}
            constants = (*cell['default_ratio'], *cell['default_angles']) if cell else (1, 1, 1, 90, 90, 90)
            fileOut.write(f'# Phase {phase}\n# MaterialName {sym.lattice or "void"}\n# Formula \n'
                          f'# Symmetry {TSL_CODES.get(sym.lattice, 0)}\n'
                          f'# LatticeConstants {" ".join(f"{i:.3f}" for i in constants)}\n')
            # ponytail: hkl families only for cubic (fcc); readers like mtex ignore them, add per lattice if OIM needs
            families = ['1 1 1', '2 0 0', '2 2 0', '3 1 1'] if sym.lattice == 'cubic' else []
            fileOut.write(f'# NumberFamilies {len(families)}\n')
            fileOut.writelines(f'# hklFamilies {hkl} 1 0.0\n' for hkl in families)
            fileOut.write('#\n')
        if ebsd.meta.get('iqScale', 1) > 1:
            fileOut.write(f'# IQ divided by {ebsd.meta["iqScale"]:g} (ebsdlab float16 range)\n#\n')
        fileOut.write(f'# GRID: {ebsd.grid}\n#\n')
        xs, ys = ebsd.xy()
        phaseID = ebsd.phaseID if ebsd.phaseID.max() > 1 else np.zeros_like(ebsd.phaseID)
        for i in range(ebsd.nPoints):
            phi1, phi, phi2 = tuple(asBungeEulers(EDAX_SAMPLE.inv() * ebsd.quaternions[i]))
            fileOut.write(f' {phi1:8.5f} {phi:8.5f} {phi2:8.5f} {xs[i]:12.5f} {ys[i]:12.5f} {ebsd.iq[i]:8.3f}'
                          f' {ebsd.ci[i]:6.3f} {phaseID[i]:2d} {int(ebsd.semSignal[i]):6d} {ebsd.fit[i]:7.3f}\n')
    print('Duration writeANG: ', int(np.round(time.time()-startTime)), 'sec')
    return

LOADERS = {'.ang': loadANG, '.osc': loadOSC, '.txt': loadTXT, '.h5': loadH5, '.crc': loadCRC, '.ctf': loadCTF}
