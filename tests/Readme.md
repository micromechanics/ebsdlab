# Help on unit testing
## For testing: execute in the home directory of the project
pytest --mpl --mpl-baseline-path=tests/baseline

## For generating the baseline images
Only do this if you are certain changes are correct

pytest --mpl-generate-path=tests/baseline tests/test_symmetry.py
