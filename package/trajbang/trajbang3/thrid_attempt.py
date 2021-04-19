#!/usr/bin/env python3

""" dans cette version on a également une accélération Ag """

import math
import cmath
import sys
import time

import sympy

check_tolerance = 1e-6

def check_is_inside(a, low, high) :
	return low - check_tolerance <= a <= high + check_tolerance

def check_is_close(a, b) :
	return abs(a - b) <= check_tolerance

def split_sign(a) :
	return abs(a), math.copysign(1.0, a)

def k(e) :
		return "\x1b[32mOK\x1b[0m" if e else "\x1b[31mKO\x1b[0m"

def c(e) :
	if 0 < e :
		return '+'
	elif e < 0 :
		return '-'
	else :
		return '0'


class TrajBang3() :

	mini_m = 0.1

	def __init__(self, jm, am, a0, s0, ag, sg) :

		self.jm, self.am, self.a0, self.s0, self.ag, self.sg = abs(jm), abs(am), a0, s0, ag, sg

		if abs(self.am) < self.mini_m :
			raise ValueError(f"Maximal acceleration can not be less than {self.mini_m}")

		if abs(self.jm) < self.mini_m :
			raise ValueError(f"Maximal jerk can not be less than {self.mini_m}")

	@property
	def val(self) :
		return { 'J_m': self.jm, 'A_m': self.am, 'A_0': self.a0, 'S_0': self.s0, 'A_g': self.ag, 'S_g': self.sg, }

	def check(self) :

		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg
		val = self.val

		def f(e) :
			try :
				return float(e.subs(val))
			except AttributeError :
				return float(e)

		cmd, dur, bch = self.compute()
		if cmd is None :
			print(f"TrajBang3({jm}, {am}, {a0}, {s0}, {ag}, {sg}) ==> KO")
			raise ValueError


		valid_condition = bch != 'Z'

		n = len(dur)
		T, J, A, S = TrajBang3.equation(n)
		for i in range(n) :
			val[f"T_{i}"] = dur[i]
			val[f"J_{i}"] = cmd[i] * val['J_m']

		# check all acceleration values are inside the limits (a0 can be out of bound)
		ai_condition = all( check_is_inside(f(i), -am, am) for i in A[1:] )

		ag_condition = check_is_close(f(A[-1]), ag)
		sg_condition = check_is_close(f(S[-1]), sg)

		total_condition = ag_condition and ai_condition and sg_condition and valid_condition

		print(f"TrajBang3({jm}, {am}, {a0}, {s0}, {ag}, {sg}) ==>", cmd, dur, bch, k( total_condition ))
		try :
			assert( total_condition )
		except AssertionError :
			print('  >', "A", [f(i) for i in A], k(ag_condition), k(ai_condition))
			print('  >', "S", [f(i) for i in S], k(sg_condition))
			raise


	@staticmethod
	def equation(n) :

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

		# print(self.val)

		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg

		# 0 step
		if sg == s0 and a0 == ag :
			return [], [], ''

		kr, mr = split_sign(ag - a0)
		qr = 0.5 * (a0 + ag) * mr / jm

		if check_is_close(sg - s0, qr) :
			return [kr], [tr], c(kr)
		else :

		


		return None, None, None


if __name__ == '__main__' :

	TrajBang3(1, 3, 0, 3, -2, 2).compute_what()



	def print_debug(u) :
		if '' in u :
			cmd, dur, res = u['']
			print((f"\x1b[32m" if res else f"\x1b[31m") + f"None\t:\x1b[0m {cmd} {dur}")

		for i in ['+0', '-0', '0+', '0-', '+-', '-+'] :
			if i in u :
				cmd, dur, res = u[i]
				print((f"\x1b[32m{i}" if res else f"\x1b[31m{i}") + f"\t:\x1b[0m {cmd} {dur}")
			else :
				print(f"\x1b[34m{i}\t: not computed\x1b[0m")

		for i in u :
			if len(i) == 3 :
				cmd, dur, res = u[i]
				print((f"\x1b[32m{i}" if res else f"\x1b[31m{i}") + f"\t:\x1b[0m {cmd} {dur}")


	u, v = TrajBang3(1, 3, 0, 3, -2, 2).compute_debug()
	print(v)
	print_debug(u)
	sys.exit(0)

	# TrajBang3(1, 2, 0, 0, 0, 6).check()
	# TrajBang3(1, 2, 0, 0, 0, -6).check()

	# # TrajBang3(1.0, 2.0, 2.0, 0.0, -1.0, -2.0).check()
	# # sys.exit(0)

	# TrajBang3(3, 3, 2, 0.0, -2, -2).check()
	# sys.exit(0)

	# TrajBang3(3.2, 3.4, 2.3, 0.3, -2.3, -1.7).check()
	
	import random

	random.seed(0)
	while True :

		jm = random.randint(5, 35) / 10.0
		am = random.randint(10, 50) / 10.0
		a0 = random.randint(-25, 25) / 10.0
		s0 = random.randint(-50, 50) / 10.0
		ag = random.randint(-25, 25) / 10.0
		sg = random.randint(-50, 50) / 10.0

		TrajBang3(jm, am, a0, s0, ag, sg).check()

		time.sleep(0.5)