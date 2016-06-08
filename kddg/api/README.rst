kddg.api
===============================================

This sub-package contains the functionality to interface with the database.

Main APIs
=========

data.py
-------

The data import API to the database. This is used to import most of the top-level data - structural and sequence information,
definition of macromolecular targets (from single chains to complexes).

db.py
-----

This in the base functionality common to all interfaces. The monomeric and protein-protein (ppi) interfaces inherit from
the classes herein.

monomer.py
----------

An interface to the monomeric stability (monomer DDG) section of the database.

Inherits from functionality in db.py.

ppi.py
----------

An interface to the binding affinity (protein-protein interface DDG) section of the database.

Inherits from functionality in db.py.

Data definition modules
=======================

schema.py
---------

This defines the relational model for the database using SQLAlchemy.


Functional modules
==================

base.py
-------

Contains basic classes.

dbi.py
------

A low-level database interface API.

A hangover from the old API. Much of this interface could be redone and some of the functionality is now in data.py.

settings.py
-----------

A module to parse and manage user-/deployment-specific settings.

High-level modules
==================

layers.py
---------

The main db.py and inherited APIs use this module to decorate functions and define a hierarchy of functionality. This module
provides the automatically-generated help functions for those APIs.
