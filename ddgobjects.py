
class Mutation(object):
	
	def __init__(self, WildTypeAA, ResidueID, MutantAA, Chain = None, SecondaryStructurePosition = None, AccessibleSurfaceArea = None):
		self.WildTypeAA = WildTypeAA
		self.ResidueID = ResidueID
		self.MutantAA = MutantAA
		self.Chain = Chain 
		self.SecondaryStructurePosition = SecondaryStructurePosition 
		self.AccessibleSurfaceArea = AccessibleSurfaceArea
	
	def __repr__(self):
		if self.Chain:
			return "%s:%s %s %s" % (self.Chain, self.WildTypeAA, self.ResidueID, self.MutantAA)
		else:
			return "?:%s %s %s" % (self.WildTypeAA, self.ResidueID, self.MutantAA)
		
	def __eq__(self, other):
		'''Only checks amino acid types and residue ID.'''
		if self.WildTypeAA != other.WildTypeAA:
			return False
		if self.ResidueID != other.ResidueID:
			return False
		if self.MutantAA != other.MutantAA:
			return False 
		return True
