from pathlib import Path
from setuptools import find_packages, setup

try:
    from cellex import __author__
except ImportError:
    __author__ = ""

setup(
    name="cellex",
    version="1.1.0",
    author=__author__,
    description="Compute single-cell cell-type expression specificity",
    long_description=Path('README.md').read_text('utf-8'),
    long_description_content_type="text/markdown",
    url="https://github.com/perslab/CELLEX",
    packages=find_packages(),
    package_data={'': ['*.txt.gz']},
    install_requires=[
        l.strip() for l in
        Path("requirements.txt").read_text("utf-8").splitlines()
    ],

    classifieres=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]


)
