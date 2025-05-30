#!/usr/bin/python3
"""TEST using the FULL set of python-requirements: create the default example that all installations create and verify it thoroughly """
import unittest
from pathlib import Path
from ebsdlab.ebsd import EBSD


class TestStringMethods(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dataDir = Path(__file__).parent/'DataFiles'

    def test_main(self):
        """
        main function
        """
        e = EBSD(str(self.dataDir/'EBSD.ang'))
        return

    def tearDown(self):
        return


if __name__ == '__main__':
    unittest.main()
