#!/usr/bin/env python3

import math

import sympy

class TrajBang3() :

	tolerance = 1e-6

	def __init__(self, jm, am, a0, s0, sg) :
		self.jm, self.am, self.a0, self.s0, self.sg = abs(jm), abs(am), a0, s0, sg

		if abs(self.am) < self.tolerance :
			raise ValueError("Maximal acceleration can not be null")

		if abs(self.jm) < self.tolerance :
			raise ValueError("Maximal acceleration can not be null")

	def check(self) :

		jm, am, a0, s0, sg = self.jm, self.am, self.a0, self.s0, self.sg

		def k(e) :
				return "\x1b[32mOK\x1b[0m" if e else "\x1b[31mKO\x1b[0m"
		def isclose(a, b) :
			return abs(a - b) < self.tolerance
		def f(e) :
			try :
				return float(e.subs(val))
			except AttributeError :
				return float(e)

		cmd, dur, bch = self.compute()

		valid_condition = bch != 'Z'

		val = { 'J_m': jm, 'A_m': am, 'A_0': a0, 'S_0': s0, 'S_g': sg, }

		n = len(dur)
		T, J, A, S = self.equation(n)
		for i in range(n) :
			val[f"T_{i}"] = dur[i]
			val[f"J_{i}"] = cmd[i] * val['J_m']

		ag_condition = isclose(A[-1].subs(val), 0.0)
		ai_condition = all( ((-am - self.tolerance) <= i.subs(val) <= (am + self.tolerance)) for i in A[1:] )
		sg_condition = isclose(S[-1].subs(val), sg)

		total_condition = ag_condition and ai_condition and sg_condition and valid_condition

		if not ( total_condition ) :
			print(f"TrajBang3({jm}, {am}, {a0}, {s0}, {sg}) ==>", cmd, dur, bch, k( total_condition ))
			print('  >', "A", [f(i) for i in A], k(ag_condition), k(ai_condition))
			print('  >', "S", [f(i) for i in S], k(sg_condition))

	def equation(self, n) :

	    T = list()
	    J = list()
	    A = [sympy.symbols('A_0'),]
	    S = [sympy.symbols('S_0'),]
	    
	    for i in range(n) :
	        T.append(sympy.symbols(f"T_{i}"))
	        J.append(sympy.symbols(f"J_{i}"))
	        A.append( (A[-1] + J[i] * T[i]).simplify() )
	        S.append( (S[-1] + A[i] * T[i] + J[i] * T[i]**2 / 2).simplify() )
	       
	    return T, J, A, S

	def compute(self) :

		jm, am, a0, s0, sg = self.jm, self.am, self.a0, self.s0, self.sg
		
		# 1 step
		k = 1.0 if 0 <= (sg - s0) else -1.0
		if (a0**2) / (2*k*jm) + s0 == sg :
			return [-k], [(a0) / (k*jm)], 'A'

		# 2 steps
		for k in [-1, 1] :
			m = (a0**2 / 2) + k*jm*(sg - s0)
			if 0 <= m :
				q = math.sqrt(m)
				d0 = (-a0 + k*q) / (k*jm)
				d1 = (q) / (jm)
				if 0 <= d0 and 0 <= d1 and q <= am :
					return [k, -k], [ d0, d1 ], 'B' + '+' if 0 < k else '-' 

		# 3 steps
		for k in [-1, 1] :
			c0 = math.copysign(1.0, k*am - a0)
			d0 = (k*am - a0) / (c0*jm)
			s01 = d0*(a0+k*am) / 2

			c2 = math.copysign(1.0, -k*am)
			d2 = (-k*am)/(c2*jm)
			s23 = d2*(k*am) / 2

			s12 = sg - s01 - s23 - s0
			d1 = s12 / (k*am)

			if 0.0 <= d1 :
				return [c0, 0, c2], [d0, d1, d2], 'C' + '+' if 0 < k else '-'

		return [0, 0, 0], [0, 0, 0], 'Z'


