#!/usr/bin/env python

"""Tests for contrasting conservation between gene clades."""

from os.path import dirname

setup_args = {}

try:
    from setuptools import setup
    # Dependencies for easy_install:
    setup_args.update(
        install_requires=[
            'Biopython >= 1.59',
            'scipy >= 0.6', 
        ])
except:
    from distutils.core import setup

setup_args.update(
    name='CladeCompare',
    version='dev',
    description=__doc__,
    author='Eric Talevich',
    author_email='etal@uga.edu',
    url='http://github.com/etal/cladecompare',
    packages=['cladecompare'],
    scripts=[(dirname(__file__) or '.') + '/cladecompare.py'],
)

setup(**setup_args)

