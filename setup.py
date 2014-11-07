#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='fr3d',
    version='0.0.1',
    packages=find_packages(include=['fr3d', 'fr3d.*']),
    url='',
    license='LICENSE.txt',
    install_requires=["numpy"],
    description='Python implementation of FR3D',
    long_description="""
    """
)
