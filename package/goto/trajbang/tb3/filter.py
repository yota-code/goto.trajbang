#!/usr/bin/env python3

class Filter() :

	def __init__(self, jm:float, am:float, debug=False) :
		self.obj = TrajBang3(jm, am)

	def reset(self, a0:float, s0:float) :
		pass

