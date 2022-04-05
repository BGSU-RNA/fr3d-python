#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='fr3d',
    version='0.0.1',
    packages=find_packages(include=['fr3d', 'fr3d.*']),
    url='',
    license='LICENSE.txt',
    install_requires=[
        'mmcif-pdbx==2.0.1',
        'numpy==1.22.3',
        'scipy==1.8.0',
        'fonttools==4.31.2',
        'matplotlib==3.5.1',
    ],
    description='Python implementation of FR3D',
    long_description="""
    """
)
