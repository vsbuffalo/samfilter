#!/usr/bin/python
from setuptools import setup
setup(
    name = "samfilter",
    version = 0.3,
    description = "Filter SAM/BAM files.",
    license = "GPL",
    scripts = ['scripts/samfilter.py'],
    install_requires = ['pysam']
    )
