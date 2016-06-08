#!/usr/bin/python2.4
# encoding: utf-8
"""
base.py
Base classes.

Created by Shane O'Connor 2012.
Copyright (c) 2012 __UCSF__. All rights reserved.
"""

import inspect

class kobject(object):

    def isOfClass(self, c):
        return self.__class__ == c

    def hasTheSameClassAs(self, o):
        return self.__class__ == o.__class__

    def getBaseClassesForObject(self):
        return self.getBaseClasses()

    def getClassName(self):
        return self.__class__.__name__

    @staticmethod
    def _getBaseClasses(c, b):
        superclasses = [sc for sc in c.__bases__ if sc not in b]
        b += superclasses
        for sc in superclasses:
            kobject._getBaseClasses(sc, b)
        return b

    @classmethod
    def getBaseClasses(cls):
        return kobject._getBaseClasses(cls, [])

    ignore_list = ["__module__", "__doc__"]

    @classmethod
    def getMembers(cls):

        lfunctions = []
        lmethods = []
        lattributes = []
        for m in inspect.getmembers(cls):
            m_name = m[0]
            m_object = m[1]
            if cls.__dict__.get(m_name):
                # Do not print inherited names
                #print(type(m_object))
                if m_name[0] != "_" and m_name not in kobject.ignore_list:
                    if inspect.isbuiltin(m_object):
                        pass
                    elif inspect.iscode(m_object):
                        pass
                    elif inspect.ismodule(m_object):
                        pass
                    elif inspect.ismethoddescriptor(m_object):
                        pass
                    elif inspect.isdatadescriptor(m_object):
                        pass
                    elif inspect.ismethod(m_object):
                        lmethods.append(m)
                    elif inspect.isfunction(m_object):
                        lfunctions.append(m)
                    elif inspect.isroutine(m_object):
                        pass
                    else:
                        lattributes.append(m)

        return {"functions" : lfunctions, "methods" : lmethods, "attributes" : lattributes}



import inspect
from klab import colortext
from string import join

class AbstractClass(object):
    '''Helper class to try to emulate abstract classes in Python. For this to work, the abstract class A must call
         super(A, self).__init__(A)
       in its __init__ method and each concrete subclass B must call
         super(B, self).__init__()
       in its __init__ method.'''

    def __init__(self, abstractClassReference, abstract_metamethods = []):
        '''Checks whether all non-private members (using the convention that private members start with a '_') are implemented.
           abstract_metamethods is an optional argument to specify non-private/metamethods which are to be considered abstract.'''
        errors = []
        for m in inspect.getmembers(abstractClassReference, inspect.ismethod):
            if m[0][0] != "_" or m[0] in abstract_metamethods:
                if getattr(self, m[0]).im_func == getattr(abstractClassReference, m[0]).im_func:
                    errors.append("The abstract method %s of class %s has not been implemented in class %s." % (m[0], abstractClassReference.__name__, self.__class__.__name__))
        if errors:
            raise NotImplementedError("\n%s" % join(errors, "\n"))


class DatasetParser(AbstractClass):
    '''I am just defining this class as an example of the proposed interface. No real OO design is happening yet.'''
    def __init__(self): super(DatasetParser, self).__init__(DatasetParser)
    def compareToDatabase(self, ProThermID, ExperimentScoreID, record): pass
    def parse(self, args = None): pass
    def addDataSet(self): pass

class AbstractDBObject(AbstractClass):
    dict = {}

    def __init__(self, ddGdb):
        self.ddGdb = ddGdb
        self.quiet = False
        self.databaseID = None
        super(AbstractDBObject, self).__init__(AbstractDBObject, ['__repr__'])

    def __getitem__(self, key):
        return self.dict[key]

    def find(self):
        '''If the records already exists in the database then two variables are returned - the record ID in the database and some printable object.
           Otherwise, false is returned.'''
        pass

    def test(self):
        DBID, details = self.find()
        if DBID != None:
            #if not self.quiet:
            #	colortext.warning("\nWarning: Experiment already exists in the database with Experiment ID=%s." % DBID)
            #	colortext.warning("*** This record ***%s\n" % self)
            #	colortext.warning("*** Database record ***\n%s\n" % details)
            return DBID

    def commit(self):
        DBID, details = self.find()
        if DBID != None:
            #if not self.quiet:
            #	colortext.error("\nExperiment already exists in the database with Experiment ID=%s." % DBID)
            #	colortext.warning("*** This record ***%s\n" % self)
            #	colortext.warning("*** Database record ***\n%s\n" % details)
            return DBID

    def remove(self):
        '''Removes the record and associated records from the database using self.databaseID.'''
        pass


class DBObject(AbstractDBObject):
    def __init__(self, ddGdb):
         super(DBObject, self).__init__(ddGdb)

    def getDatabaseID(self):
        if not self.databaseID:
            raise Exception("Cannot get the database ID of an uncommitted record.")
        else:
            return self.databaseID

