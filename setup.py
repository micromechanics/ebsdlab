#!/usr/bin/env python
""" setup to generate version number"""
from __future__ import annotations

from setuptools import setup

import releaseVersion

if __name__ == '__main__':
    setup(name='ebsdlab', version=releaseVersion.getVersion()[1:])
