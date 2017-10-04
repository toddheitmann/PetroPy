"""Setup script for PetroPy"""

from setuptools import setup
from os import path
from petropy import __version__

with open(path.join(path.dirname(__file__), "requirements.txt"), "r") as f:
    requirements = f.read().splitlines()

setup(
    name = 'petropy',
    packages=["petropy", ],
    version = __version__,
    description = 'A package to calculate petrophysical properties for formation evaluation.',
    author = 'Todd Heitmann',
    author_email = 'toddheitmann@protonmail.com',
    url = 'https://github.com/toddheitmann/petropy',
    keywords = ['petrophysics', 'formation evaluation', 'reservoir characterization'],
    classifiers=[
        "Intended Audience :: Customer Service",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Other Audience",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: System :: Filesystems",
        "Topic :: Scientific/Engineering :: Information Analysis",
    ],
    install_requires = requirements
)
