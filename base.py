import inspect

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
	