#!/usr/bin/env python

"""Test for sites of contrasting conservation between gene clades."""

from os.path import dirname

setup_args = {}

try:
    from setuptools import setup
    # Dependencies for easy_install:
    setup_args.update(
        install_requires=[
            'biofrills >= 0.2',
            'biopython >= 1.60',
            'scipy >= 0.6',
            'reportlab >= 2.5',
            'weblogo >= 3.0',
            'bottle >= 0.10',
        ])
except ImportError:
    from distutils.core import setup


DIR = (dirname(__file__) or '.') + '/'


setup_args.update(
    name='CladeCompare',
    version='0.2',
    description=__doc__,
    author='Eric Talevich',
    author_email='etal@uga.edu',
    url='http://github.com/etal/cladecompare',
    packages=['cladecomparelib'],
    scripts=[
        DIR + 'cladecompare.py',
        DIR + 'cladereport.py',
        DIR + 'cladeweb.py'
    ])

setup(**setup_args)

