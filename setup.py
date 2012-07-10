#!/usr/bin/python
from setuptools import setup
import samfilter

setup(
    name = "samfilter",
    version = samfilter.__version__,
    description = "Filter SAM/BAM files.",
    license = "GPL",
    scripts = ['scripts/samfilter.py'],
    install_requires = ['pysam']
    )
