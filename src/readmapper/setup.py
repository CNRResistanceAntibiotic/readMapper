# setup.py file install libraries wanted for readmapper
# launch it with pip install -e .
from setuptools import setup, find_packages

setup(name='readmapper', version='1.0', packages=find_packages(), install_requires=['pandas', 'docx'])
