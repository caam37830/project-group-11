# setup file
from setuptools import setup
import sys, os
import setuptools

this_dir = os.path.dirname(os.path.realpath(__file__))

__version__ = '0.1.1'


setup(
    name='sir',
    version=__version__,
    author='James Keane, Song Liang, Luke Chen',
    author_email='jmkeane@uchicago.edu',
    description='package with tools for analysing SIR model',
    url="https://github.com/caam37830/project-group-11",
    python_requires='>=3.6',
    packages=['sir'],
    zip_safe=True,

)