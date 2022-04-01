#!/usr/bin/env python3

import goto.trajbang.tb3.compute

class Tb3Filter() :

	def __init__(self, jm:float, am:float, a0=0.0, s0=0.0, period=0.1) :
		self.obj = goto.trajbang.tb3.compute.Tb3Compute(jm, am)
		
		self.a0, self.s0 = a0, s0
		self.period = period

	def reset(self, a0:float, s0:float) :
		pass

	def run(self, sg, ag) :
		self.ag, self.sg = ag, sg

		j = self._jerk()

		a1 = self.a0 + self.period * j
		s1 = self.s0 + self.period * (self.a0 + a1) / 2

		self.a0, self.s0 = a1, s1

		return s1, a1, j

	def _jerk(self) :
		# try to compute a direct solution
		j0 = (( -(self.ag + (3.0 * self.a0)) * self.period ) + 2.0 * (self.sg - self.s0)) / (2.0 * self.period**2)
		j1 = (( ((3.0 * self.ag) + self.a0) * self.period ) - 2.0 * (self.sg - self.s0)) / (2.0 * self.period**2)

		a1 = self.a0 + j0 * self.period

		if ( -self.obj.jm <= j0 <= self.obj.jm ) and ( -self.obj.jm <= j1 <= self.obj.jm ) and ( -self.obj.am <= a1 <= self.obj.am ) :
			# if the direct solution is compatible with the limits
			return j0

		r_lst = self.obj._compute(self.a0, self.s0, self.ag, self.sg)

		t = 0.0
		c = 0.0
		for cmd, dur in r_lst :
			if (t + dur) <= self.period :
				c += cmd * dur
			else :
				c += cmd * max(self.period - t, min(0.0, self.period))
			t += dur

		return c * self.obj.jm / self.period

if __name__ == '__main__' :
	import numpy as np
	import matplotlib.pyplot as plt

	x = [0.0 for i in range(8)] + [10.0 for i in range(700)]
	u = Tb3Filter(0.01, 0.25)
	y = np.array([u.run(i, 0.0) for i in x])

	plt.plot(x)
	plt.plot(y[:,0])
	plt.show()