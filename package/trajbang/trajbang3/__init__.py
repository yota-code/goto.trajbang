#!/usr/bin/env python3

import math

import sympy

class TrajBang3() :

	def __init__(self, jm, am, a0, s0, sg) :
		print(f"TrajBang3({jm}, {am}, {a0}, {s0}, {sg})")
		self.jm, self.am, self.a0, self.s0, self.sg = jm, am, a0, s0, sg

	def check(self) :

		def ok(e) :
				return "\x1b[32mOK\x1b[0m" if e else "\x1b[31mKO\x1b[0m"
		def isclose(a, b) :
			return math.isclose(a, b, rel_tol=1e-5, abs_tol=1e-5)
		def f(e) :
			try :
				return float(e.subs(val))
			except AttributeError :
				return float(e)

		cmd, dur, bch = self.compute()
		val = {
			'J_m': self.jm,
			'A_m': self.am, 'A_0': self.a0,
			'S_0': self.s0, 'S_g': self.sg,
		}

		n = len(dur)
		T, J, A, S = self.equation(n)
		for i in range(n) :
			val[f"T_{i}"] = dur[i]
			val[f"J_{i}"] = cmd[i] * val['J_m']

		print(' --> ', cmd, dur, bch,
			"A", [f(i) for i in A], ok( isclose(A[-1].subs(val), 0.0) ),
			"S", [f(i) for i in S], ok( isclose( S[-1].subs(val), self.sg) )
		)

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
		ws = 1.0 if 0 <= (sg - s0) else -1.0
		if (-a0**2) / (-2*jm*ws) + s0 == sg :
			return [-ws], [(-a0) / (jm*k)], "A"

		# 2 steps
		# try with first solution
		k = -1
		dn = (a0**2/2) + k*jm*(sg - s0)
		if 0 <= dn :
			qn = math.sqrt(dn)
			if qn <= am and -qn <= a0 :
				return [k, -k], [ (-a0 - qn) / (jm * k), (-qn) / (jm * k) ], "B-"

		# try with second solution
		k = 1
		dp = (a0**2/2) + k*jm*(sg - s0)
		if 0 <= dp :
			qp = math.sqrt(dp)
			if qp <= am and a0 <= qp :
				return [k, -k], [ (-a0 + qp) / (jm * k), (qp) / (jm * k) ], "B+"

if __name__ == '__main__' :
	TrajBang3(1, 2, 2, 0, -1).check()
	TrajBang3(1, 2, -2, 0, 1).check()
