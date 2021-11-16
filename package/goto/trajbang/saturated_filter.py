#!/usr/bin/env python3

import collections

period = 0.001

def bound(x, b) :
	return max(-b, min(x, b))

class Integral1() :
	def __init__(self, x1=0.0, y1=0.0) :
		self.x1 = x1
		self.y1 = y1

	def step(self, x0) :
		self.x0 = x0
		self.y0 = self.y1 + period * (self.x0 + self.x1) / 2.0

		self.x1 = self.x0
		self.y1 = self.y0

		return self.y0

class Integral1_sat() :
	def __init__(self, sat, x1=0.0, y1=0.0) :
		self.x1 = x1
		self.y1 = y1

		self.sat = abs(sat)

	def step(self, x0) :
		self.x0 = x0
		self.y0 = bound( self.y1 + period * (self.x0 + self.x1) / 2.0, self.sat )

		self.x1 = self.x0
		self.y1 = self.y0

		return self.y0

class Filter2() :
	def __init__(self, epsilon, omega0) :
		self.epsilon = epsilon
		self.omega0 = omega0

		self.xref1_prev, self.xref0_prev = 0.0, 0.0

		self.integral1_1 = Integral1()
		self.integral1_2 = Integral1()

		self.hist = collections.defaultdict(list)

	def step(self, xcom) :
		self.xcom = xcom

		self.xref2 = self.omega0**2 * ( xcom - self.xref0_prev - (self.xref1_prev * 2 * self.epsilon / self.omega0) )
		self.xref1 = self.integral1_2.step(self.xref2)
		self.xref0 = self.integral1_1.step(self.xref1)

		self.hist['xref2'].append(self.xref2)
		self.hist['xref1'].append(self.xref1)
		self.hist['xref0'].append(self.xref0)

		self.xref0_prev = self.xref0
		self.xref1_prev = self.xref1

		return self.xref0

class Filter2_sat() :
	def __init__(self, epsilon, omega0, sat1, sat2) :

		self.epsilon = epsilon
		self.omega0 = omega0

		self.xref1_prev, self.xref0_prev = 0.0, 0.0

		self.integral1_1 = Integral1()
		self.integral1_2 = Integral1_sat(sat1)

		self.sat1 = abs(sat1)
		self.sat2 = abs(sat2)

		self.hist = collections.defaultdict(list)

	def step(self, xcom) :
		self.xcom = xcom

		self.xref2 = bound( self.omega0**2 * ( xcom - self.xref0_prev - (self.xref1_prev * 2 * self.epsilon / self.omega0) ), self.sat2 )
		self.xref1 = self.integral1_2.step(self.xref2)
		self.xref0 = self.integral1_1.step(self.xref1)

		self.hist['xref2'].append(self.xref2)
		self.hist['xref1'].append(self.xref1)
		self.hist['xref0'].append(self.xref0)

		self.xref0_prev = self.xref0
		self.xref1_prev = self.xref1

		return self.xref0
