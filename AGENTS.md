# Agent Notes

## Project Structure
- `ebsdlab/` contains the package code for EBSD file import, orientation/quaternion math, symmetry handling, and plotting.
- `tests/` contains pytest tests and image-comparison baselines. Test data is under `tests/DataFiles/`.
- `docs/` contains the Sphinx documentation.
- Packaging is configured with `setup.cfg` and `setup.py`; dependencies are pinned in `requirements.txt` and `requirements-dev.txt`.

## Environment
- The project supports Python >=3.10.
- A local `.venv/` may exist in the repository. Use it when it is ready, but do not recreate, delete, or modify it without checking with the user.
- Development dependencies are listed in `requirements-dev.txt`.

## Commands
- Install locally: `python -m pip install .`
- Install development dependencies when approved: `python -m pip install -r requirements-dev.txt`
- Run tests from the repository root: `pytest --mpl --mpl-baseline-path=tests/baseline`
- Build docs: `make -C docs html`
- Type check: `python -m mypy ebsdlab`
- Lint: `python -m pylint ebsdlab`

## Conventions
- Keep scientific behavior changes small and covered by focused tests.
- Do not regenerate baseline images unless the visual change is intentional and reviewed.
- Preserve user data files and local environment directories.
- Prefer NumPy/SciPy APIs over ad hoc numerical code when they make behavior clearer.
- Keep imported OSC Euler angles and scalar data as `float16` unless a task explicitly
  requires higher precision. This is a deliberate memory/performance trade-off for
  large maps: EBSD indexing angular accuracy is typically no better than about 0.1°,
  so retaining `float32` input precision does not usually improve the physical result.

## Issues
- Open issues are tracked in GitHub Issues for `micromechanics/ebsdlab`.
- Local issue notes may also appear in repository documentation when GitHub is not the only working context.
