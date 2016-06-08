kddg - The Kortemme Lab DDG Database Repository
===============================================

Database APIs for a database used by our lab to store data about mutageneses, macromolecular complexes, experimental measurements, and computational methods. This database is used primarily for DDG predictions at present.

Coding convention
=================

- General rules      : The PEP 8 conventions (https://www.python.org/dev/peps/pep-0008/) are mostly sensible.
- Blank lines        : See "Blank Lines" in PEP 8. Particularly, surround top-level function and class definitions with two blank lines.
- Line length        : Can be more than 79 characters long. We have widescreen technology.
- Whitespace         : Use spaces instead of tabs, with indents set to 4 spaces.
- Natural language   : Please do not use contractions.
- Nerdy references   : Allowed / encouraged.
- Shiny?             : Shiny! See above.

Installation
============

Once the PyPI package is up and running (it is currently not)...

This package can be installed via:
::
  pip install kddg

To install via pip and allow git push/pulling, use:

In a virtualenv:
::
  pip install -e git+ssh://git@github.com/Kortemme-Lab/kddg.git#egg=Package

In your user-directory:
::
  pip install --user -e git+ssh://git@github.com/Kortemme-Lab/kddg.git#egg=Package

In your user-directory, without using SSH:
::
  pip install --user -e git+https://github.com/Kortemme-Lab/kddg.git#egg=Package


