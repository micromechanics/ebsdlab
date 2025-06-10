#!/usr/bin/env python
""" setup to generate version number"""
from __future__ import annotations
import configparser
from setuptools import setup

if __name__ == '__main__':
    config = configparser.ConfigParser()
    config.read('setup.cfg')
    setup(name='ebsdlab', version=config['metadata']['version'])
