#!/usr/bin/env python2

import distutils.core
import subprocess, shlex

# Uploading to PyPI
# =================
# The first time only:
# $ python setup.py register -r pypi
#
# Every version bump:
# $ git tag <version>; git push --tags
# $ python setup.py sdist upload -r pypi

version = '0.29'
distutils.core.setup(
    name='kddg',
    version=version,
    author='Kortemme Lab, UCSF',
    author_email='support@kortemmelab.ucsf.edu',
    url='https://github.com/Kortemme-Lab/kddg',
    download_url='https://github.com/Kortemme-Lab/kddg/tarball/'+version,
    license='MIT',
    description="Database APIs for managing mutageneses, macromolecular complexes, experimental measurements, and computational methods for our lab.",
    long_description=open('README.rst').read(),
    keywords=['database', 'api', 'biophysics'],
    packages=['kddg'],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 3 - Alpha",
    ],
)
