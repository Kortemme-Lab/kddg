#!/usr/bin/python2.4
# encoding: utf-8
"""
score.py

THIS FILE SHOULD BE REMOVED ONCE THE DDG SCHEDULER HAS BEEN REWRITTEN.

Created by Shane O'Connor 2012.
Copyright (c) 2012 __UCSF__. All rights reserved.
"""

class ddgScore(object):
	type = None
	version = None
	
	def __init__(self, data = {}):
		self.data = data
	
	def getType(self):
		return type
	
	def getVersion(self):
		return version
	
	def setData(self, data):
		self.data = data

	def getScores(self):
		return {"type" : self.type, "version" : self.version, "data" : self.data}

class ddgTestScore(ddgScore):
	type = "test"
	version = "0.1"