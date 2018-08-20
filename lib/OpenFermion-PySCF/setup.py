from setuptools import setup, find_packages

# This reads the __version__ variable from openfermionpyscf/_version.py
exec(open('openfermionpyscf/_version.py', encoding='utf-8').read())

# Readme file as long_description:
long_description = open('README.rst', encoding='utf-8').read()

# Read in requirements.txt
requirements = open('requirements.txt', encoding='utf-8').readlines()
requirements = [r.strip() for r in requirements]

setup(
    name='openfermionpyscf',
    version=__version__,
    author='The OpenFermion Developers',
    author_email='help@openfermion.org',
    url='http://www.openfermion.org',
    description='A plugin allowing OpenFermion to interface with PySCF.',
    long_description=long_description,
    install_requires=requirements,
    license='Apache 2',
    packages=find_packages()
)
