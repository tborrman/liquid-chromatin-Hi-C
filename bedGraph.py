import numpy as np

class bedGraph:
	'''
	Class for maniuplating bedGraphs
	'''

	def __init__(self, bG_file):
		self.name = bG_file[:-9]
		self.FH = open(bG_file, 'r')
		self.chr = []
		self.start = []
		self.end = []
		self.values = np.array([])

		for line in self.FH:
			splitline=line.strip().split('\t')
			self.chr.append(splitline[0])
			self.start.append(int(splitline[1]))
			self.end.append(int(splitline[2]))
			self.values = np.append(self.values, float(splitline[3]))






	



