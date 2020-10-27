#! /usr/bin/env python3 

import setuptools
import os
import re

PACKAGE = setuptools.find_packages(exclude=('tests',))[0]

def get_readme():
    with open("README.md", "r") as fh:
        text = fh.read()
    return text

def read_init(splitlines=False):
    here = os.path.abspath(os.path.dirname(__file__))
    filename = os.path.join(here, PACKAGE, '__init__.py')
    with open(filename) as fh:
        lines = fh.read()
    if splitlines:
        lines = lines.splitlines()
    return lines

def get_version():
    lines = read_init(splitlines=True)
    for line in lines:
        match = re.search('__version__\s*=\s*([\'"])(.*)\\1', line)
        if match:
            return match.groups()[1]
    raise RuntimeError("Unable to find version string.")

# Python package install
 
setuptools.setup(
    name='pyoifits',
    version=get_version(),
    packages=setuptools.find_packages(),
    license='LICENSE.txt',
    author="Régis Lachaume",
    author_email="regis.lachaume@gmail.com",
    description='OIFITS (Data Exchange Standard for Optical Interferometry, A&A 597, A8, 2017)',
    long_description=get_readme(),
    long_description_content_type='text/markdown',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: Public Domain",
    ],
    python_requires='>=3.8',
    install_requires=[
        "astropy>=4.0", 
        "numpy>=1.17",
        "scipy>=1.5",
        "astroquery>=0.4",
    ],
)
