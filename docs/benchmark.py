"""Resource comparison of ebsdlab, orix and MTEX for Python for docs/source/comparisonLibraries.rst.

Each package runs in its own virtual environment and every measurement in a fresh process, so imports and peak
memory do not mix. Each file runs twice per package, sharing an empty numba cache (NUMBA_CACHE_DIR): where numba
compiles, the first run includes the compilation, as after an installation, and the times read "first / second".
Usage:

    python docs/benchmark.py EBSDLAB_VENV ORIX_VENV MTEX_VENV FILE [FILE ...]
    python docs/benchmark.py --plot-only

EBSDLAB_VENV holds only `pip install .`, ORIX_VENV only `pip install orix`, MTEX_VENV only
`pip install "mtex[numba,plot,io]"`. A measurement writes docs/benchmark_results.json and plots the first FILE to
docs/source/_static/benchmark_results.png, which docs/source/comparisonLibraries.rst shows; --plot-only redraws the
figure from the existing JSON. The data files are not part of the repository.
"""
import json
import os
import platform
import subprocess
import sys
from datetime import date
from pathlib import Path
from tempfile import TemporaryDirectory

RESULTS = Path(__file__).parent/'benchmark_results.json'
FIGURE = Path(__file__).parent/'source'/'_static'/'benchmark_results.png'


def worker(package: str, fileName: str) -> None:
    """Load one file and plot its IPF map along Z; print the measurements as JSON."""
    import resource  # pylint: disable=import-outside-toplevel
    import time  # pylint: disable=import-outside-toplevel
    from importlib.metadata import version  # pylint: disable=import-outside-toplevel
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
    print(json.dumps({'version': version(package), 'python': platform.python_version(), 'points': int(points),
                      'load': timeLoad, 'plot': timePlot, 'rssImport': rssImport,
                      'rssPeak': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                      'data': nbytes(data)/1024**2}))


def main() -> None:
    """Run the workers in all environments, print the results, save them as JSON and plot them."""
    venvs = {'ebsdlab': Path(sys.argv[1]), 'orix': Path(sys.argv[2]), 'mtex': Path(sys.argv[3])}
    results: dict = {'date': date.today().isoformat(), 'venvGB': {}, 'files': []}
    for package, venv in venvs.items():
        size = sum(f.stat().st_size for f in venv.rglob('*') if f.is_file() and not f.is_symlink())
        results['venvGB'][package] = size/1024**3
        print(f'{package} venv: {size/1024**3:.2f} GB')
    for fileName in sys.argv[4:]:
        sizeMB = Path(fileName).stat().st_size/1024**2
        print(f'\n{Path(fileName).name}, {sizeMB:.0f} MB on disk')
        packages = {}
        for package, venv in venvs.items():
            with TemporaryDirectory() as cache:
                runs = [json.loads(subprocess.run([venv/'bin'/'python', __file__, '--worker', package, fileName],
                                                  capture_output=True, text=True, check=True,
                                                  env=os.environ | {'NUMBA_CACHE_DIR': cache}
                                                  ).stdout.strip().splitlines()[-1]) for _ in range(2)]
                compiled = any(Path(cache).iterdir())
            packages[package] = {'compiled': compiled, 'runs': runs}
            result = runs[-1]

            def times(key: str) -> str:
                """Time of the second run, preceded by the first one where numba compiled."""
                if result[key] is None:
                    return 'not possible'
                return (f'{runs[0][key]:.1f} / ' if compiled else '') + f'{result[key]:.1f} s'
            print(f"  {package:8s} {result['points']:>9,d} points  load {times('load')}  "
                  f"IPF map {times('plot')}  peak memory {result['rssPeak']:6.0f} MB "
                  f"(after imports {result['rssImport']:4.0f} MB)  loaded data {result['data']:6.1f} MB")
        results['files'].append({'name': Path(fileName).name, 'sizeMB': sizeMB, 'packages': packages})
    RESULTS.write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
    plot(results)


def plot(results: dict) -> None:
    """Bar charts of the first file, one panel per measurement: smaller is better.

    Where numba compiled, the solid bar is the later run and the light bar extends it to the first run.
    """
    import matplotlib  # pylint: disable=import-outside-toplevel
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel
    packages = results['files'][0]['packages']
    names = {'ebsdlab': 'ebsdlab', 'orix': 'orix', 'mtex': 'MTEX'}
    colors = {'ebsdlab': '#2a78d6', 'orix': '#eb6834', 'mtex': '#1baf7a'}
    panels = [('load', 'time to load the data', 's', '.1f'), ('plot', 'time to plot the IPF map', 's', '.1f'),
              ('rssPeak', 'peak RAM', 'MB', '.0f'), ('rssImport', 'RAM after importing the package', 'MB', '.0f'),
              ('data', 'RAM of the loaded map', 'MB', '.0f'), ('venvGB', 'installed size on disk', 'GB', '.2f')]
    fig, axes = plt.subplots(2, 3, figsize=(12, 4), layout='constrained')
    fig.get_layout_engine().set(hspace=0.12)
    for iAx, (ax, (key, title, unit, fmt)) in enumerate(zip(axes.flat, panels)):
        values = {}
        for package, entry in packages.items():
            later = results['venvGB'][package] if key == 'venvGB' else entry['runs'][-1][key]
            first = entry['runs'][0][key] if key in ('load', 'plot') and entry['compiled'] else None
            values[package] = (first, later)
        right = max(max(v for v in pair if v is not None) for pair in values.values() if pair[1] is not None)
        for i, (package, (first, later)) in enumerate(values.items()):
            if later is None:
                ax.text(0, i, ' not possible', va='center', color='0.4')
                continue
            ax.barh(i, later, height=0.9, color=colors[package])
            label = f'{later:{fmt}} {unit}'
            if first is not None:
                if first > later:
                    ax.barh(i, first-later, left=later, height=0.9, color=colors[package], alpha=0.35)
                label = f'{first:{fmt}} / {later:{fmt}} {unit}'
            ax.text(max(later, first or 0)+0.02*right, i, label, va='center')
        ax.set_yticks(range(len(values)), [names[p]+('*' if packages[p]['compiled'] else '') for p in values] if iAx % 3 == 0 else [])
        if iAx % 3 == 0:
            ax.get_yticklabels()[list(values).index('ebsdlab')].set_fontweight('bold')
        ax.invert_yaxis()
        ax.set_xlim(0, 1.45*right)
        ax.set_xticks([])
        ax.tick_params(left=False)
        ax.spines[['top', 'right', 'bottom']].set_visible(False)
        ax.set_title(title, loc='left')
    fig.suptitle('Resource usage: smaller is better', fontsize='xx-large')
    fig.supxlabel('* light bar: first run with numba compilation', x=0.01, ha='left', fontsize='medium')
    fig.savefig(FIGURE, dpi=150)


if __name__ == '__main__':
    if sys.argv[1] == '--worker':
        worker(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == '--plot-only':
        plot(json.loads(RESULTS.read_text(encoding='utf-8')))
    else:
        main()
