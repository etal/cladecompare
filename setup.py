#!/usr/bin/env python

"""Tests for contrasting conservation between gene clades."""

try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(  name='CladeCompare',
        version='dev',
        description=__doc__,
        author='Eric Talevich',
        author_email='eric.talevich@gmail.com',
        url='http://github.com/etal/cladecompare',
        packages=['cladecompare'],
        scripts=['cladecompare.py'],
        )

