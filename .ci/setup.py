from setuptools import setup, find_packages

setup(
    name='pdbx',
    version='0.0.1',
    author='NDB',
    author_email='',
    packages=find_packages(include=['pdbx', 'pdbx.*']),
    url='',
    license='LICENSE.txt',
    description='Pure python tools to work with cif data from PDB',
    long_description='Pure python tools to work with cif data from PDB',
)
