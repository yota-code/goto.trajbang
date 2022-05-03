#!/usr/bin/env python3

import random
import sys

from trajbang.trajbang3 import TrajBang3

for i in range(int(sys.argv[1])) :
	print(i, end='\r')
	vec = [round(random.uniform(-128.0, 128.0), 1) for i in range(5)]
	if vec[0] == 0.0 or vec[1] == 0.0 :
		continue
	if i % 100 == 0 :
		print(vec)
	TrajBang3(* vec).check()