#!/usr/bin/env python
import argparse
import bedGraph as bG
import numpy as np
import copy

parser=argparse.ArgumentParser(description='Input bedGraph A and bedGraph B and output bedGraph of log2(A/B)')
parser.add_argument('-a', help='input bedGraph A', type=str, required=True)
parser.add_argument('-b', help='input bedGraph B', type=str, required=True)
args=parser.parse_args()

def check_coordinates(x, y):
	'''
	Given bedGraph objects x and y
	test if coordinates are matching
	Args:
	x : bedGraph object
	y : bedGraph object
	Returns:
	b = boolean True if coordinates are equal
	'''
	b = (x.chr == y.chr) and (x.start == y.start) and (x.end == y.end)
	return b

def log2ratio(x,y):
	'''
	Given bedGraph objects x and y
	return bedGraph of log2ratio(x,y)
	Args:
	x : bedGraph object
	y : bedGraph object
	Returns:
	z : bedGraph object of log2ratio(x/y)
	'''
	lr = np.log2(x.values/y.values)
	z = copy.deepcopy(x)
	z.values = lr
	return z

	

def main():

	a = bG.bedGraph(args.a)
	b = bG.bedGraph(args.b)
	if check_coordinates(a,b) :
		c = log2ratio(a,b)
		for i in range(len(c.start)):
			print '\t'.join([c.chr[i], str(c.start[i]), str(c.end[i]), str(round(c.values[i], 4))])

if __name__ == '__main__':
	main()


