from distutils.core import setup

setup(
    name='pdbx',
    version='0.0.1',
    author='NDB',
    author_email='',
    packages=['pdbx', 'pdbx.alternate', 'pdbx.bird', 'pdbx.cc', 'pdbx.core',
              'pdbx.ndb', 'pdbx.pdb', 'pdbx.persist', 'pdbx.reader',
              'pdbx.writer'],
    url='',
    license='LICENSE.txt',
    description='Pure python tools to work with cif data from PDB',
    long_description=open('README.mkd').read(),
)
