#!/usr/bin/env python
from setuptools import setup

setup(name="DynPy",
      version="0.1",
      description="Module for engineering calculations on dynamic systems",
      packages=["dynpy"],
      install_requires=["pylatex", "matplotlib", "timer", "pypandoc"])
