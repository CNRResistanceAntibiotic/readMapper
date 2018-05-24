#!/usr/bin/python3

# setup.py file install libraries wanted for readmapper
# launch it with pip install -e .
from setuptools import setup, find_packages
import glob

setup(
    name='readmapper',
    version='1.0.1',
    description='ReadMapper: pipeline CNR Resistance with Ariba tool',
    packages=find_packages(),
    author='Richard Bonnet',
    author_email='rbonnet@chu-clermontferrand.fr',
    url='https://github.com/CNRResistanceAntibiotic/readMapper',
    scripts=glob.glob('scripts/*'),
    install_requires=['pandas', 'docx'],
    license='GPLv3',
    classifiers=[
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],

)
