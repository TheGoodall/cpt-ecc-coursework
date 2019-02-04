import math
import numpy as np

from hammingGeneratorMatrix import *

def message(a):
	r = 2
	while 2**(r) -2*r -1 < len(a):
		r += 1
	k = 2**(r) -r -1
	print("k", k)
	print("r", r)
	out = []
	
	out += (decimalToVector(len(a), r))
	out += a
	while len(out) < k:
		out.append(0)
	return out	


