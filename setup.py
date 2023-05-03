from setuptools import setup, find_packages

from leechem import __version__

setup(
    name='leechem',
    version=__version__,

    url='https://github.com/wagenadl/leechem',
    author='Daniel A. Wagenaar',
    author_email='daw@caltech.edu',
    packages=find_packages(),
    install_requires=['numpy']
    #py_modules=['leechem'],
)
