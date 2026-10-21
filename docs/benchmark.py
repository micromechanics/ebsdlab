"""Resource comparison of ebsdlab, orix and MTEX for Python for docs/source/comparisonOrix.rst.

Each package runs in its own virtual environment and every measurement in a fresh process, so imports and peak
memory do not mix. Usage:

    python docs/benchmark.py EBSDLAB_VENV ORIX_VENV MTEX_VENV FILE [FILE ...]

EBSDLAB_VENV holds only `pip install .`, ORIX_VENV only `pip install orix`, MTEX_VENV only
`pip install "mtex[numba,plot,io]"`. A one-time comparison: the results,
with versions and date, are in docs/source/comparisonOrix.rst; the data files are not part of the repository.
"""
import json
import subprocess
import sys
from pathlib import Path


def worker(package: str, fileName: str) -> None:
    """Load one file and plot its IPF map along Z; print the measurements as JSON."""
    import resource  # pylint: disable=import-outside-toplevel
    import time  # pylint: disable=import-outside-toplevel
    import matplotlib  # pylint: disable=import-outside-toplevel
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel
    import numpy as np  # pylint: disable=import-outside-toplevel

    def nbytes(obj, seen=None, depth=0) -> int:
        """Bytes of all numpy arrays reachable from obj (attributes, dicts, lists), each array counted once."""
        seen = set() if seen is None else seen
        if id(obj) in seen or depth > 5:
            return 0
        seen.add(id(obj))
        if isinstance(obj, np.ndarray):
            return obj.nbytes
        if type(obj).__name__ == 'Rotation' and hasattr(obj, 'as_quat'):  # scipy Rotation
            return obj.as_quat().nbytes
        children = list(obj.values()) if isinstance(obj, dict) else list(obj) if isinstance(obj, (list, tuple)) else []
        children += list(getattr(obj, '__dict__', {}).values())
        children += [getattr(obj, n, None) for c in type(obj).__mro__ for n in getattr(c, '__slots__', ())
                     if isinstance(n, str)]
        return sum(nbytes(child, seen, depth+1) for child in children)

    if package == 'ebsdlab':
        from ebsdlab.ebsd import EBSD  # pylint: disable=import-outside-toplevel
    elif package == 'mtex':
        import mtex  # pylint: disable=import-outside-toplevel
    else:
        from orix import io  # pylint: disable=import-outside-toplevel
        from orix.plot import IPFColorKeyTSL  # pylint: disable=import-outside-toplevel
        from orix.vector import Vector3d  # pylint: disable=import-outside-toplevel
    rssImport = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
    start = time.perf_counter()
    if package == 'ebsdlab':
        data = EBSD(fileName)
    elif package == 'mtex':
        data = mtex.loadEBSD(fileName)
    else:
        data = io.load(fileName)
    timeLoad = time.perf_counter()-start
    start = time.perf_counter()
    timePlot: float | None = None
    try:
        if package == 'ebsdlab':
            data.plotIPF('ND', show=False)
        elif package == 'mtex':
            rgb = np.full(data.shape+(3,), np.nan)
            for phaseID in data.indexedPhaseId:
                phase = data[data.mineralList[phaseID]]
                key = mtex.ipfTSLKey(phase.CS, ipfDirection=mtex.vector3d.Z)
                rgb[data.phaseId == phaseID] = key.orientation2color(phase.orientations)
            mtex.plot(data, rgb)
        else:
            rgb = np.zeros((data.size, 3))
            for phaseID, phase in data.phases_in_data:
                key = IPFColorKeyTSL(phase.point_group, direction=Vector3d.zvector())
                rgb[data.phase_id == phaseID] = key.orientation2color(data[phase.name].orientations)
            data.plot(rgb, return_figure=True)
        plt.gcf().canvas.draw()
        timePlot = time.perf_counter()-start
    except IndexError:  # orix maps need a square grid; hex grids fail in CrystalMap.get_map_data
        pass
    points = data.nPoints if package == 'ebsdlab' else (~data.isPadding).sum() if package == 'mtex' else data.size
    print(json.dumps({'points': int(points),
                      'load': timeLoad, 'plot': timePlot, 'rssImport': rssImport,
                      'rssPeak': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                      'data': nbytes(data)/1024**2}))


def main() -> None:
    """Run the workers in all environments and print the results."""
    venvs = {'ebsdlab': Path(sys.argv[1]), 'orix': Path(sys.argv[2]), 'mtex': Path(sys.argv[3])}
    for package, venv in venvs.items():
        size = sum(f.stat().st_size for f in venv.rglob('*') if f.is_file() and not f.is_symlink())
        print(f'{package} venv: {size/1024**3:.2f} GB')
    for fileName in sys.argv[4:]:
        print(f'\n{Path(fileName).name}, {Path(fileName).stat().st_size/1024**2:.0f} MB on disk')
        for package, venv in venvs.items():
            output = subprocess.run([venv/'bin'/'python', __file__, '--worker', package, fileName],
                                    capture_output=True, text=True, check=True).stdout
            result = json.loads(output.strip().splitlines()[-1])
            plot = 'not possible' if result['plot'] is None else f"{result['plot']:5.1f} s"
            print(f"  {package:8s} {result['points']:>9,d} points  load {result['load']:6.1f} s  "
                  f"IPF map {plot}  peak memory {result['rssPeak']:6.0f} MB "
                  f"(after imports {result['rssImport']:4.0f} MB)  loaded data {result['data']:6.1f} MB")


if __name__ == '__main__':
    if sys.argv[1] == '--worker':
        worker(sys.argv[2], sys.argv[3])
    else:
        main()
