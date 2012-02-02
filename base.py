class kobject(object):
	
	def isOfClass(self, c):
		return self.__class__ == c

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