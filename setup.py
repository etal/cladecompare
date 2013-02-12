#!/usr/bin/env python

"""Test for sites of contrasting conservation between gene clades."""

from os.path import dirname

setup_args = {}

try:
    from setuptools import setup
    # Dependencies for easy_install:
    setup_args.update(
        install_requires=[
            'biofrills >= 0.1',
            'biopython >= 1.58',
            'scipy >= 0.6', 
            'reportlab >= 2.5', 
        ])
except:
    from distutils.core import setup


DIR = (dirname(__file__) or '.') + '/'


setup_args.update(
    name='CladeCompare',
    version='dev',
    description=__doc__,
    author='Eric Talevich',
    author_email='etal@uga.edu',
    url='http://github.com/etal/cladecompare',
    packages=['cladecompare'],
    scripts=[DIR + 'cladecompare.py', DIR + 'cladereport.py'],
)

setup(**setup_args)

