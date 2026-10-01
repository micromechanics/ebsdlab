# Agent Notes

## Project Structure

- `ebsdlab/` contains the package code for EBSD file import, orientation/quaternion math, symmetry handling, and plotting.
- `tests/` contains pytest tests and image-comparison baselines. Test data is under `tests/DataFiles/`.
- `docs/` contains the Sphinx documentation.
- Packaging is configured with `setup.cfg` and `setup.py`; dependencies are pinned in `requirements.txt` and `requirements-dev.txt`.

## Commands

- Install locally: `python -m pip install .`
- Install development dependencies when approved: `python -m pip install -r requirements-dev.txt`
- Run tests: `pytest --mpl --mpl-baseline-path=tests/baseline`
- Build docs: `make -C docs html`
- Type check: `python -m mypy ebsdlab`
- Lint: `python -m pylint ebsdlab`

## Conventions

- The global rules in `~/.codex/AGENTS.md` apply, with these local decisions:
  - Two empty lines above every `def`, also inside classes (except `__init__`).
  - Single quotes for strings; docstrings keep triple double quotes.
  - Keep useful existing comments; move them with their code when reordering. Doxygen markers (`##`, `# @file`,
    `# @brief`, `# @name`, `# @{`, `# @}`) and commented-out code may be removed.
  - Leave tests as they are; add or change tests only when asked. `pytest --mpl` is a visual check: run it only
    when asked.

- Do not regenerate baseline images unless the visual change is intentional and reviewed.
- Baseline review, when a change alters images:
  1. `rm -rf tests/baseline_review && pytest --mpl --mpl-baseline-path=tests/baseline --mpl-results-path=tests/baseline_review --mpl-generate-summary=html`
     writes `baseline.png`, `result.png` and `result-failed-diff.png` per failing test, plus `fig_comparison.html`.
  2. Add a side-by-side `tests/baseline_review/side_by_side.png` (matplotlib, baseline left, result right, one row
     per changed test); for new behaviour not covered by a test, add panels showing it, with timings if speed changed.
  3. Look at the images yourself, then report what differs and why; do not commit `tests/baseline_review/`.
  4. Only after approval copy each `result.png` to `tests/baseline/<test>.png` and to `docs/source/_static/` where a
     copy exists there; rerun the tests.
- Prefer NumPy/SciPy APIs over ad hoc numerical code when they make behavior clearer.
- `phaseID` is `uint8`: 0 = not identified (`self.sym[0]` is `Symmetry()`), phases are numbered from 1 and
  `self.sym[k]` is the symmetry of phase `k`. Loaders convert file conventions to this; `writeANG` converts back.
- Keep imported OSC Euler angles and scalar data as `float16` unless a task explicitly requires higher precision. This is a deliberate memory/performance trade-off for large maps: EBSD indexing angular accuracy is typically no better than about 0.1°, so retaining `float32` input precision does not usually improve the physical result.

## Issues

- Open issues are tracked in GitHub Issues for `micromechanics/ebsdlab`.
- Local issue notes may also appear in repository documentation when GitHub is not the only working context.
