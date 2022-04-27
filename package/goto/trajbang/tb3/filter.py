#!/usr/bin/env python3

import goto.trajbang.tb3.compute

class Tb3Filter() :

	def __init__(self, jm:float, am:float, a0=0.0, s0=0.0, period=0.1) :
		self.obj = goto.trajbang.tb3.compute.Tb3Compute(jm, am)
		
		self.a0, self.s0 = a0, s0
		self.period = period

	def reset(self, a0:float, s0:float) :
		pass

	def run(self, ag, sg) :
		self.ag, self.sg = ag, sg

		j = self._jerk()

		a1 = self.a0 + self.period * j
		s1 = self.s0 + self.period * (self.a0 + a1) / 2

		self.a0, self.s0 = a1, s1

		return s1, a1, j

	def _jerk(self) :
		# j0 = (( -(self.ag + (3.0 * self.a0)) * self.period ) + 2.0 * (self.sg - self.s0)) / (2.0 * self.period**2)
		epsilon = 1e-6

		# a 1 step solution
		j0 = (self.ag - self.a0) / self.period
		s1 = self.period * (self.ag + self.a0) / 2.0 + self.s0

		print("=== ", self.a0, self.s0, self.ag, self.sg)

		print("quick ?", j0, s1, self.sg)
		if abs(s1 - self.sg) < epsilon :
			print("quick !")
			return j0

		# a 2 step solution
		j0 = (( -(self.ag + (3.0 * self.a0)) * self.period ) + 2.0 * (self.sg - self.s0)) / (2.0 * self.period**2)
		j1 = (( ((3.0 * self.ag) + self.a0) * self.period ) - 2.0 * (self.sg - self.s0)) / (2.0 * self.period**2)
		a1 = j0 * self.period + self.a0

		print(f"!!! a={self.a0} -> {self.ag}\n    s={self.s0} -> {self.sg}\n    {j0}, {j1}, {a1}")

		if ( -self.obj.jm - epsilon <= j0 <= self.obj.jm + epsilon) and ( -self.obj.jm - epsilon <= j1 <= self.obj.jm + epsilon ) and ( -self.obj.am - epsilon <= a1 <= self.obj.am + epsilon ) :
			print("express", j0, j1, a0)
			return j0

		r_lst = self.obj._compute(self.a0, self.s0, self.ag, self.sg)
		print("normal", r_lst)
		# 	# if the direct solution is compatible with the limits
		# 	print("return quick:", j0)
		# 	return j0

		t = 0.0
		c = 0.0
		for cmd, dur in r_lst :
			if (t + dur) <= self.period :
				c += cmd * dur
			else :
				c += cmd * max(self.period - t, min(0.0, self.period))
			t += dur
		j = c * self.obj.jm / self.period

		return j

if __name__ == '__main__' :
	import numpy as np
	import matplotlib.pyplot as plt

	from cc_pathlib import Path

	a0, s0, ag, sg = 0.0, 70.0, 0.0, 55.0
	x_lst = [s0 for i in range(8)] + [sg for i in range(200)]

	u = Tb3Filter(1.0/3.0, 4.0/3.0, 0.0, 70.0)
	u_comp = u.obj._compute(a0, s0, ag, sg)
	u_inte = u.obj._integrate(a0, s0, ag, sg, u_comp)

	v = Tb3Filter(1.0/3.0, 4.0/3.0, 0.0, 70.0)
	v_comp = u.obj._compute(a0, s0, ag, sg)
	v_inte = u.obj._integrate(a0, s0, ag, sg, v_comp)

	y_arr = np.array([u.run(i, 0.0) for i in x_lst])
	z_lst = list()
	for n, i in enumerate(y_arr) :
		print("--> ", n)
		z_lst.append(v.run(i[0], i[1]))
	z_arr = np.array(z_lst)

	plt.figure()
	plt.subplot(3,1,1)
	plt.plot(y_arr[:,0])
	plt.plot(z_arr[:,0])
	plt.subplot(3,1,2)
	plt.plot(y_arr[:,1])
	plt.plot(z_arr[:,1])
	plt.subplot(3,1,3)
	plt.plot(y_arr[:,2])
	plt.plot(z_arr[:,2])
	plt.show()

	print(y_arr)
	print(z_arr)

	m_lst = [
		["#000_ref", "001_dot", "002_ref", "003_dot", "004_is_requested", "005_value", "006_mode", "007_is_limited", "008__I3_spd_fin", "009__I4_distance_remaining"],
	]
	for t, (s, a, j) in enumerate(y_arr) :
		m_lst.append([s, a, 0.0, 0.0, 1 if (50 <= t < 70) else 0, 75.0, 1, 0, 55.0, u_inte['D'][-1][0]])

	Path("input.tsv").save(m_lst)
		