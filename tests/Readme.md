# Help on unit testing
## For testing: execute in the home directory of the project
pytest --mpl --mpl-baseline-path=tests/baseline

## For generating the baseline images
Only do this if you are certain changes are correct

pytest --mpl-generate-path=tests/baseline tests/test_symmetry.py

## Reviewing image differences
Review failing images before regenerating baselines; many failures are only 1-2 px layout shifts.

```
pytest --mpl --mpl-baseline-path=tests/baseline --mpl-results-path=/tmp/mpl-review --mpl-generate-summary=html
python tests/shift_review.py /tmp/mpl-review /tmp/shift-review [maxShift=2] [tolerance=0.1]
python -m http.server 8765 --directory /tmp   # open http://127.0.0.1:8765/shift-review/index.html
```

`shift_review.py` marks a changed pixel green if a consistent shift of at most `maxShift` px makes its
7x7 neighbourhood match the baseline, otherwise red. It shows a baseline/result blink view, the overlay,
a shift map and a x5 zoom on the reddest area (before / overlay / after). Nothing is approved automatically.

- Green only: layout shift, data unchanged.
- Red outlines of text: sub-pixel text rendering, harmless.
- Red ring on circles or axes: scaling; rerun with a larger `maxShift` (e.g. 6).
- Red inside the data area: real change; fix the code, not the baseline.

In the pytest-mpl diff image black means identical; moved text and lines appear white.

Regenerate only the approved tests, then `git add` new baseline files:
```
pytest --mpl-generate-path=tests/baseline tests/test_ebsd.py::test_ebsd_ci
```
