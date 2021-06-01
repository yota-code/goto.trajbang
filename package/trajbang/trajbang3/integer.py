#!/usr/bin/env python3

""" dans cette version on a également une accélération Ag """

import cmath
import collections
import fractions
import itertools
import math
import random
import sys
import time

import sympy
import numpy as np
import matplotlib.pyplot as plt

from cc_pathlib import Path

trajbang_T = collections.namedtuple('trajbang_T', ['cmd', 'dur'])

check_tolerance = 1e-6

def check_is_inside(a, low, high) :
	return low - check_tolerance <= a <= high + check_tolerance

def check_is_close(a, b) :
	return abs(a - b) <= check_tolerance

def split_value_sign(a) :
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


def egcd(a, b):
    if a == 0:
        return b, 0, 1
    else:
        gcd, x, y = egcd(b % a, a)
        return gcd, y - (b // a) * x, x

print( egcd(30, 50) )


class TrajBang3() :

	mini_m = 0.1

	def __init__(self, jm: float, am: float, a0: float, s0: float, ag: float, sg: float, tm: float) :

		self.tm = fractions.Fraction(int(1000 * tm), 1000) # real period of the computation in seconds
		self.jm = fractions.Fraction(int(1000 * abs(jm)), 1000)
		self.am = fractions.Fraction(int(1000 * abs(am)), 1000)

		print(self.tm, self.jm, self.am)

		self.a0, self.s0 = a0, s0
		self.ag = min(self.am, max(- self.am, ag))
		self.sg = sg

		self.compute_kj()

	def compute_kj(self) :
		jm, tm = self.jm, self.tm
		print(jm, tm, jm * tm, 6 * (jm * tm).denominator)

		v = float(jm * tm)

		for i in range(1, 3000) :
			if (v * i) % 6 == 0 :
				print(i)
				break


		self.kt = 1.0
		self.kj = 1.0

	@property
	def val(self) :
		return { 'J_m': self.jm, 'A_m': self.am, 'A_0': self.a0, 'S_0': self.s0, 'A_g': self.ag, 'S_g': self.sg, }

	def integrate(self, result_lst, plot=None) :
		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg

		T = [0,]
		J = [0,]
		A = [a0,]
		S = [s0,]
		D = [0,]

		for cmd, dur in result_lst :
			t = int(dur)

			T.append(T[-1] + t)
			J.append(int(cmd) * jm)

			A.append( (A[-1] + J[-1]*t) )
			S.append( (S[-1] + A[-2]*t + J[-1]*t**2 / 2) )
			D.append( (D[-1] + S[-2]*t + A[-2]*t**2 / 2 + J[-1]*t**3 / 6) )

		print(f"J = {', '.join(str(i) for i in J)}")
		print(f"A = {', '.join(str(i) for i in A)}")
		print(f"S = {', '.join(str(i) for i in S)}")
		print(f"D = {', '.join(str(i) for i in D)}")

		if plot :
			plt.subplot(4, 1, 1)
			plt.step(T, J)
			plt.subplot(4, 1, 2)
			plt.plot(T, A)
			plt.subplot(4, 1, 3)
			plt.plot(T, S)
			plt.subplot(4, 1, 4)
			plt.plot(T, D)
			plt.show()

	def check(self) :

		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg
		val = self.val

		def f(e) :
			try :
				return float(e.subs(val))
			except AttributeError :
				return float(e)

		res = self.compute()

		cmd = [cmd for cmd, dur in res]
		dur = [dur for cmd, dur in res]

		n = len(dur)

		T = list()
		J = list()
		A = [a0,]
		S = [s0,]
		D = [0.0,]

		for i in range( len(dur) ) :
			T.append(dur[i])
			J.append(cmd[i] * jm)
			A.append( (A[-1] + J[i] * T[i]) )
			S.append( (S[-1] + A[i] * T[i] + J[i] * T[i]**2 / 2) )
			D.append( (D[-1] + S[i] * T[i] + A[i] * T[i]**2 / 2 + J[i] * T[i]**3 / 6) )

		print(D)

		# check all acceleration values are inside the limits (a0 can be out of bound)
		ai_condition = all( check_is_inside(f(i), -am, am) for i in A[1:] )

		ag_condition = check_is_close(f(A[-1]), ag)
		sg_condition = check_is_close(f(S[-1]), sg)

		total_condition = ag_condition and ai_condition and sg_condition

		print(f"TrajBang3(jm={jm}, am={am}, a0={a0}, s0={s0}, ag={ag}, sg={sg}) ==>", k( total_condition ))

		print("> cmd :", '\t'.join(f"{i:5g}" for i in cmd))
		print("> dur :", '\t'.join(f"{i:5g}" for i in dur))
		print("-" * 4)
		print("  > T :", '\t'.join(f"{i:5g}" for i in ([0.0,] + list(itertools.accumulate(T)))))
		print("  > A :", '\t'.join(f"{i:5g}" for i in A), f" \t=( {ag} )>", k(ag_condition), k(ai_condition))
		print("  > S :", '\t'.join(f"{i:5g}" for i in S), f" \t=( {sg} )>", k(sg_condition))
		#	raise

		if not total_condition :
			with Path('error_listing.log').open('at') as fid :
				fid.write(f"TrajBang3(jm={jm}, am={am}, a0={a0}, s0={s0}, ag={ag}, sg={sg}) ==>\n\t{cmd}\n\t{dur}\n")

	def get_q(self, a_from, a_to) :
		m, w = split_value_sign(a_to - a_from)
		t = m / self.jm
		q = (a_to + a_from) * t / 2
		return q, t, w		

	def compute(self) :

		res = np.zeros((8, 2))

		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg

		# computation of the initial segment, from A0 to the closest of Am or -Am

		ap = min(max(a0, -am), am)
		# if am < a0:
		# 	ap = am
		# elif a0 < -am :
		# 	qi, ti, wi = self.get_q(a0, -am)
		# 	ap = -am
		# else :
		# 	qi, ti, wi = 0.0, 0.0, 0.0
		# 	ap = a0
		qi, ti, wi = self.get_q(a0, ap)

		res[0] = [wi, ti]
		print(f"qi={qi} ti={ti} wi={wi} ap={ap}")

		#computation of the final segment, from Ap to Ag
		qr, tr, wr = self.get_q(ap, ag)

		res[4] = [wr, tr]
		print(f"qr={qr} tr={tr} wr={wr}")

		# computation of the remaining part to complete up to qd
		qd = sg - s0 - qi - qr
		print(f"qd={qd}")

		ab = ag if ( 0 <= qd * wr ) else ap

		print(f"ab={ab}")

		# if qd != 0.0 :
		wd = math.copysign(1.0, qd)

		#print(f"wd={wd}")

		de = ab**2 + jm*abs(qd)
		ad_1 = ( -ab + math.sqrt(de) ) * wd
		ad_2 = ( -ab - math.sqrt(de) ) * wd

		if 0 <= ad_1 :
			ad = ad_1
		else :
			ad = ad_2

		#ab = abs(an)

		print(f"ad_1={ad_1} ad_2={ad_2} ad={ad} wd={wd}")

		i = 4 if ( 0 <= qd * wr ) else 0

		at = ab + ad*wd

		print(f"at={at}")

		if am < abs(at) :
			res[i+1] = [wd, abs(am*wd - ab)/jm]
			res[i+2] = [0.0, (at**2 - am**2)/(am*jm)]
			res[i+3] = [-wd, abs(am*wd - ab)/jm]
		else :
			res[i+1] = [wd, ad / jm]
			res[i+3] = [-wd, ad / jm]

			# res[i+1] = [wd, ad / jm]
			# res[i+3] = [-wd, ad / jm]

		#else :
		#	ad, wd = 0.0, 0.0

		return res


if __name__ == '__main__' :

	result_lst = [
		[6, 20],
		[0, 20],
		[-6, 20]
	]
	u = TrajBang3(1.4, 3.7, 0, 0, 0, 8, 0.08)
	sys.exit(0)
	u.integrate(result_lst, True)
	sys.exit(0)

	# TrajBang3(jm=20, am=8, a0=53, s0=-93, ag=4, sg=-24).check()
	# TrajBang3(jm=1, am=2, a0=3, s0=0, ag=1, sg=2.75).check()
	# TrajBang3(jm=1, am=2, a0=3, s0=0, ag=1, sg=-0.5).check()
	# TrajBang3(jm=30, am=45, a0=40, s0=-15, ag=5, sg=-35).check()
	TrajBang3(1, 2, 0, 0, 0, 8).check()
	# TrajBang3(1, 2, 0, 0, 0, 3).check()
	sys.exit(0)

	TrajBang3(1, 2, 3, 0, 1, -0.5).check()
	# TrajBang3(jm=6, am=9, a0=8, s0=-3, ag=1, sg=-7).check()

	
	test_vec = [
		[1, 2, 0, 0, 0, 9],
		[1, 2, 1, -0.5, 0, 0],
		[1, 1, 0.5, 0, 0, 0.5],
		[1, 1, 0.5, 0, 0, 0.5],
		[1, 1, 0.5, 0, 0, 0.5],
		[1, 1, 0.5, 0, 0, 0.5],
		[1, 1, 0, 0, 0, 0.5],
		[1, 1, 0, 0, 0, 0.5],
		[1, 1, 0, 0, 0, 2],
		[1, 1, 0, 0, 0, -2],
		[1, 1, -0.25, 0, 0, -2],
		[1, 1, 0.25, 0, 0, 2],
		[1, 1, 0.75, 0, 0, -2],
		[1, 1, -0.75, 0, 0, 2],
		[1, 2.0, -1.0, -2.0, 0, 4.0],
		[1, 2, 0, 0, 0, 6],
		[1, 2, 3, -2, 0, 3],
		[1, 2, 3, -2, 0, 6],
		[1, 2, 3, 0, -1, 2],
		[1, 2, 0, 0, 0, 10],
		[1, 2, 0, 0, 0, 9],
		[1, 2, 3, 0, 1, 2.75],
		[1, 2, 3, 0, 1, -0.5]
	]

	for vec in test_vec :
		TrajBang3(* vec).check()
	sys.exit(0)

	m = 128.0
	while True:	
		jm = round(random.uniform(0.1, m), 1)
		am = round(random.uniform(0.1, m), 1)
		a0 = round(random.uniform(-m, m), 1)
		s0 = round(random.uniform(-m, m), 1)
		ag = round(random.uniform(-am, am), 1)
		sg = round(random.uniform(-m, m), 1)

		TrajBang3(jm, am, a0, s0, ag, sg).check()

	sys.exit()



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