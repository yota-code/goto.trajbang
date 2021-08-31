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
		
		# 1 step
		t0, k0 = split_sign(ag - a0)
		s1 = a0*t0 + k0*jm*t0**2 / 2 + s0
		if check_is_close(s1, sg) :
			return [k0], [t0], c(k0)

		# 3 steps
		for k in [-1, 1] :
			t0, k0 = split_sign((k*am - a0) / jm)
			s01 = a0*t0 + k0*jm*t0**2 / 2

			t2, k2 = split_sign((ag - k*am) / jm)
			s23 = k*am*t2 + k2*jm*t2**2 / 2

			t1 = (sg - s01 - s23 - s0) / (k*am)

			# print(f"k={k}, s01={s01}, s23={s23}, t={[t0, t1, t2]}")

			if 0 <= t1 :
				return [k0, 0, k2], [t0, t1, t2], ''.join(c(i) for i in [k0, k, k2])

		# 2 steps
		ad = (ag - a0) / jm
		aq = (ag**2 - a0**2)
		sd = jm * (sg - s0)

		if ag != 0 :
			t0_pz = ad
			t1_pz = (-aq + 2*sd) / (2*ag*jm)
			if 0 < t0_pz and 0 < t1_pz :
				return [1, 0], [t0_pz, t1_pz], '+0'

			t0_mz = -ad
			t1_mz = (aq + 2*sd) / (2*ag*jm)
			if 0 < t0_mz and 0 < t1_mz :
				return [-1, 0], [t0_mz, t1_mz], '-0'
		else :
			t0_pz, t1_pz, t0_mz, t1_mz = None, None, None, None

		if a0 != 0 :
			t0_zp = (-aq + 2*sd) / (2*a0*jm)
			t1_zp = ad
			if 0 < t0_zp and 0 < t1_zp :
				return [0, 1], [t0_zp, t1_zp], '0+'

			t0_zm = (aq + 2*sd) / (2*a0*jm)
			t1_zm = -ad
			if 0 < t0_zm and 0 < t1_zm :
				return [0, 1], [t0_zm, t1_zm], '0-'
		else :
			t0_zp, t1_zp, t0_zm, t1_zm = None, None, None, None

		ap = (ag**2 + a0**2)

		qp = ap/2 + sd
		if 0 <= qp :
			t0_pm = -(a0 + math.sqrt(qp)) / jm
			t1_pm = -(ag + math.sqrt(qp)) / jm
			if 0 < t0_pm and 0 < t1_pm :
				return [1, -1], [t0_pm, t1_pm], '+-'
		else : 
			t0_pm, t1_pm = None, None

		qd = ap/2 - sd
		if 0 <= qd :
			t0_mp = (a0 - math.sqrt(qd)) / jm
			t1_mp = (ag - math.sqrt(qd)) / jm
			if 0 < t0_mp and 0 < t1_mp :
				return [-1, 1], [t0_mp, t1_mp], '-+'
		else : 
			t0_mp, t1_mp = None, None

		return None, None, None


	def compute_debug(self) :

		result_map = dict()

		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg

		# 0 step

		result_map[''] = ( [], [], sg == s0 and a0 == ag )
		
		# 1 step
		t0, k0 = split_sign(ag - a0)
		s1 = a0*t0 + k0*jm*t0**2 / 2 + s0
		result_map[c(k0)] = ( [k0], [t0], check_is_close(s1, sg) )

		# 3 steps
		for k in [-1, 1] :
			t0, k0 = split_sign((k*am - a0) / jm)
			s01 = a0*t0 + k0*jm*t0**2 / 2

			t2, k2 = split_sign((ag - k*am) / jm)
			s23 = k*am*t2 + k2*jm*t2**2 / 2

			t1 = (sg - s01 - s23 - s0) / (k*am)

			print(''.join(c(i) for i in [k0, k, k2]), t1)

			result_map[''.join(c(i) for i in [k0, k, k2])] = ( [k0, 0, k2], [t0, t1, t2], 0 <= t1 )

		# 2 steps
		ad = (ag - a0) / jm
		aq = (ag**2 - a0**2)
		sd = jm * (sg - s0)

		if ag != 0 :
			t0_pz = ad
			t1_pz = (-aq + 2*sd) / (2*ag*jm)
			result_map['+0'] = ( [1, 0], [t0_pz, t1_pz], 0 < t0_pz and 0 < t1_pz )

			t0_mz = -ad
			t1_mz = (aq + 2*sd) / (2*ag*jm)

			result_map['-0'] = ( [-1, 0], [t0_mz, t1_mz], 0 < t0_mz and 0 < t1_mz )

		if a0 != 0 :
			t0_zp = (-aq + 2*sd) / (2*a0*jm)
			t1_zp = ad
			result_map['0+'] = ( [0, 1], [t0_zp, t1_zp], 0 < t0_zp and 0 < t1_zp )

			t0_zm = (aq + 2*sd) / (2*a0*jm)
			t1_zm = -ad
			result_map['0-'] = ( [0, 1], [t0_zm, t1_zm], 0 < t0_zm and 0 < t1_zm )

		ap = (ag**2 + a0**2)

		qp = ap/2 + sd
		if 0 <= qp :
			t0_pm = -(a0 + math.sqrt(qp)) / jm
			t1_pm = -(ag + math.sqrt(qp)) / jm
			result_map['+-'] = ( [1, -1], [t0_pm, t1_pm], 0 < t0_pm and 0 < t1_pm )

		qd = ap/2 - sd
		if 0 <= qd :
			t0_mp = (a0 - math.sqrt(qd)) / jm
			t1_mp = (ag - math.sqrt(qd)) / jm
			result_map['-+'] = ( [-1, 1], [t0_mp, t1_mp], 0 < t0_mp and 0 < t1_mp )

		return result_map, self.val

	def compute_what(self) :
		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg

		k2, m2 = split_sign(ag - a0)

		# the time to transit from a0 to ag with a variation of ±jm
		tr = m2 / jm
		sr = a0 * tr + 0.5 * k2 * jm * tr**2

		k1, m1 = split_sign(sg - s0 - sr)



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