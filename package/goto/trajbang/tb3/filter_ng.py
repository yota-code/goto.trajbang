#!/usr/bin/env python3


class Intergral_0() :
	def __init__(self, period=0.05) :
		self.period = period
		self.x_prev = 0.0

	def run(self, x) :
		y = self.period * self.x_prev / (1 - self.x_prev)
		self.x_prev = x