#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

test_requirements = ['pytest>=3', ]
setup_requirements = ['pytest-runner', ]
install_requires = ['numpy', 'matplotlib', 'scipy', 'scikit-learn', 'pytest']

setup(
    name='ebsdlab',
    author='Steffen Brinckmann',
    author_email='sbrinckm@gmail.com',
    python_requires='>=3',
    classifiers=[
        'Development Status :: 3 - Beta',
        'Intended Audience :: Experts',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    description='Read, process and plot EBSD (electron backscatter diffraction) data',
    install_requires=['numpy', 'matplotlib', 'scipy'],
    license='GNU General Public License v3',
    long_description=readme,
    include_package_data=True,
    keywords='EBSD',
    packages=find_packages('ebsdlab', exclude=['*tests*']),
    package_dir={'': 'ebsdlab'},
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/micromechanics/ebsdlab.git',
    version='0.0.2',
    zip_safe=False,
)
