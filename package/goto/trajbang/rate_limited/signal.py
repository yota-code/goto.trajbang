#!/usr/bin/env python3

dt = 10 # ms

class RateLimiter() :
	def __init__(self, upper_limit, lower_limit=None) :
		self.upper_limit = upper_limit
		self.lower_limit = -upper_limit if lower_limit is None else lower_limit

		self.prev = 0.0

	def run(self, x) :
		d = max(self.lower_limit * dt / 1000.0, min(x - self.prev, self.upper_limit * dt / 1000.0))
		y = self.prev + d
		self.prev = y
		return y

class Integrate() :
	def __init__(self) :
		self.prev = 0.0

	def run(self, x) :
		y = (x * dt / 1000.0) + self.prev
		self.prev = y
		return y

class Derivative() :
	def __init__(self) :
		self.prev = 0.0

	def run(self, x) :
		y = (x - self.prev) / (dt / 1000.0)
		self.prev = x
		return y

class JerkLimiter() :
	def __init__(self, am, jm) :
		self.jm = jm
		self.am = am

		self.jrl = RateLimiter(self.jm)
		self.arl = RateLimiter(self.am)
		self.deriv = Derivative()
		self.integ = Integrate()

	def run(self, x) :
		y1 = self.arl.run(x)
		y2 = self.deriv.run(y1)
		y3 = self.jrl.run(y2)
		y4 = self.integ.run(y3)
		return y4

