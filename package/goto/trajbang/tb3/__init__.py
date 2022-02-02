#!/usr/bin/env python3

import itertools
import math
import random
import sys
import time

import numpy as np
import matplotlib.pyplot as plt
import sympy

from cc_pathlib import Path

def split_value_sign(a) :
	return abs(a), math.copysign(1.0, a) if a != 0.0 else 0.0

def sign_of(a) :
	if a > 0 :
		return 1
	elif a < 0 :
		return -1
	else :
		return 0
	

def sym_eval(e, v) :
	try :
		return float(e.subs(v))
	except AttributeError :
		return float(e)


class Trajectory() :

	mini_m = 0.1
	check_tolerance = 1e-6

	def __init__(self, jm: float, am: float, debug=False) :
		self.jm, self.am = abs(jm), abs(am)

		if abs(self.am) < self.mini_m :
			raise ValueError(f"Maximal acceleration can not be less than {self.mini_m}")

		if abs(self.jm) < self.mini_m :
			raise ValueError(f"Maximal jerk can not be less than {self.mini_m}")

		self.debug = debug

	def _prepare(self, a0: float, s0: float, ag: float, sg: float) :
		self.a0, self.s0, self.ag, self.sg = a0, s0, max(-self.am, min(ag, self.am)), sg
		self.res = np.zeros((8, 2))

		return self.jm, self.am, self.a0, self.s0, self.ag, self.sg

	def check_is_inside(self, a, low, high) :
		return low - self.check_tolerance <= a <= high + self.check_tolerance

	def check_is_close(self, a, b) :
		return abs(a - b) <= self.check_tolerance

	@property
	def val(self) :
		return { 'J_m': self.jm, 'A_m': self.am, 'A_0': self.a0, 'S_0': self.s0, 'A_g': self.ag, 'S_g': self.sg, }

	@staticmethod
	def _integrate_sym(length=3) :
		""" symbol integration """
		Jm, Am, A0, S0, T = sympy.symbols('J_m A_m A_0 S_0 T')

		T = sympy.symbols(' '.join(f'T_{i}' for i in range(length)))
		J = sympy.symbols(' '.join(f'J_{i}' for i in range(length)))

		A = [A0,]
		S = [S0,]
		D = [0,]
		for i in range(length) :
			A.append( (A[-1] + J[i] * T[i]) )
			S.append( (S[-1] + A[i] * T[i] + J[i] * T[i]**2 / 2) )
			D.append( (D[-1] + S[i] * T[i] + A[i] * T[i]**2 / 2 + J[i] * T[i]**3 / 6) )

		return A, S, D

	def _integrate_num(self) :
		""" numerical integration """

		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg
		t0, d0 = 0.0, 0.0

		self.poly = { k : list() for k in ['T', 'A', 'S', 'D'] }

		self.poly['T'].append(t0)
		
		for i, (cmd, dur) in enumerate(self.sol) :
			t = dur
			self.poly['T'].append(self.poly['T'][-1] + t)

			a1 = cmd * jm
			self.poly['A'].append([a1, a0])

			s2 = a1 / 2.0
			s1 = a0
			self.poly['S'].append([s2, s1, s0])

			d3 = s2 / 3.0
			d2 = s1 / 2.0
			d1 = s0
			self.poly['D'].append([d3, d2, d1, d0])

			a0 = a1 * t + a0
			s0 = s2 * t**2 + s1 * t + s0
			d0 = d3 * t**3 + d2 * t**2 + d1 * t + d0

		self.poly['A'].append([a0,])
		self.poly['S'].append([s0,])
		self.poly['D'].append([d0,])

		return self.poly

	def get_at_time(self, t) :
		""" return the A, S and D reached at the time passed in argument,
		outside the bounds, boundary values are returned
		 """
		t_start = self.poly['T'][0]
		t_stop = self.poly['T'][-1]
		# print(f'start[{t_start}] <= {t} < stop[{t_stop}]')
		if t < t_start :
			return (self.a0, self.s0, 0.0)
		if t_stop <= t :
			return (self.poly['A'][-1][-1], self.poly['S'][-1][-1], self.poly['D'][-1][-1])
		for i in range(len(self.sol)) :
			t_before = self.poly['T'][i]
			t_after = self.poly['T'][i+1]
			# print(f'  {i}. before[{t_before}] < {t} <= after[{t_after}]')
			if t_before <= t < t_after :
				t = t - t_before
				a1, a0 = self.poly['A'][i]
				s2, s1, s0 = self.poly['S'][i]
				d3, d2, d1, d0 = self.poly['D'][i]
				return (
					a1,
					a1 * t + a0,
					s2 * t**2 + s1 * t + s0,
					d3 * t**3 + d2 * t**2 + d1 * t + d0
				)

	@property
	def title(self) :
		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg
		return f"TrajBang3(jm={jm}, am={am}, a0={a0}, s0={s0}, ag={ag}, sg={sg})"

	def plot(self, result_dir=None) :
		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg
		cmd, dur, T, J, A, S, D = self.cmd, self.dur, self.T, self.J, self.A, self.S, self.D

		plt.figure()
		plt.subplot(4, 1, 1)
		plt.step(T, J)
		plt.subplot(4, 1, 2)
		plt.plot(T, A)
		plt.subplot(4, 1, 3)
		plt.plot(T, S)
		plt.subplot(4, 1, 4)
		plt.plot(T, D)
		if result_dir is None :
			plt.show()
		else :
			plt.savefig(result_dir / f"{self.title}.png")		

	def status(self) :
		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg
		cmd, dur, T, J, A, S, D = self.cmd, self.dur, self.T, self.J, self.A, self.S, self.D
		print("cmd :", ' '.join(f"{i:7.3g}" for i in cmd))
		print("dur :", ' '.join(f"{i:7.3g}" for i in dur))
		print("-" * 4)
		print("  T :", ' '.join(f"{i:7.3g}" for i in ([0.0,] + list(itertools.accumulate(T)))))
		print("  A :", ' '.join(f"{i:7.3g}" for i in A + [ag,]))
		print("  S :", ' '.join(f"{i:7.3g}" for i in S + [sg,]))
		print("  D :", ' '.join(f"{i:7.3g}" for i in D))

	def check(self, a0: float, s0: float, ag: float, sg: float, result_dir=None) :
		jm, am, a0, s0, ag, sg = self._prepare(a0, s0, ag, sg)

		val = self.val

		def f(e) :
			try :
				return float(e.subs(val))
			except AttributeError :
				return float(e)

		res = self.compute(a0, s0, ag, sg)

		cmd, dur, T, J, A, S, D = self._integrate_sym(res)

		# check all acceleration values are inside the limits (except a0 which can be out of bounds)
		ai_condition = all(self.check_is_inside(f(i), -am, am) for i in A[1:] )

		ag_condition = self.check_is_close(f(A[-1]), ag)
		sg_condition = self.check_is_close(f(S[-1]), sg)

		total_condition = ag_condition and ai_condition and sg_condition

		print(f"TrajBang3(jm={jm}, am={am}, a0={a0}, s0={s0}, ag={ag}, sg={sg}) ==>", k( total_condition ))

		print("> cmd :", '\t'.join(f"{i:5g}" for i in cmd))
		print("> dur :", '\t'.join(f"{i:5g}" for i in dur))
		print("-" * 4)
		print("  > T :", '\t'.join(f"{i:5g}" for i in ([0.0,] + list(itertools.accumulate(T)))))
		print("  > A :", '\t'.join(f"{i:5g}" for i in A), f" \t=( {ag} )>", k(ag_condition), k(ai_condition))
		print("  > S :", '\t'.join(f"{i:5g}" for i in S), f" \t=( {sg} )>", k(sg_condition))
		print("  > D :", '\t'.join(f"{i:5g}" for i in D))

		if not total_condition :
			with Path('error_listing.log').open('at') as fid :
				fid.write(f"\n{self.title} ==>\n\t{cmd}\n\t{dur}\n")

	def get_q(self, a_from, a_to) :
		m, w = split_value_sign(a_to - a_from)
		t = m / self.jm
		q = (a_to + a_from) * t / 2
		return q, t, w

	def compute(self, a0: float, s0: float, ag: float, sg: float) :

		jm, am, a0, s0, ag, sg = self._prepare(a0, s0, ag, sg)

		# computation of the initial segment, from A0 to Ap the closest of either Am or -Am
		ap = max(-am, min(a0, am))
		qi, ti, wi = self.get_q(a0, ap)

		self.res[0] = [wi, ti]
		if self.debug :
			print(f"ap={ap}")
			print(f"qi={qi} ti={ti} wi={wi}")

		#computation of the final segment, from Ap to Ag
		qr, tr, wr = self.get_q(ap, ag)

		self.res[4] = [wr, tr]
		if self.debug :
			print(f"qr={qr} tr={tr} wr={wr}")

		# computation of Qd the remaining part to go from S0 to Sg without the initial and final segments
		qd = sg - s0 - qi - qr
		if self.debug :
			print(f"qd={qd}")

		ab = ag if ( 0 <= qd * wr ) else ap
		if self.debug :
			print(f"ab={ab}")

		# if Qd != 0.0 we need to compute the remaining command
		wd = math.copysign(1.0, qd)

		de = ab**2 + jm*abs(qd)
		ad_1 = ( -ab + math.sqrt(de) ) * wd
		ad_2 = ( -ab - math.sqrt(de) ) * wd
		ad = ad_1 if 0 <= ad_1 else ad_2
		if self.debug :
			print(f"ad_1={ad_1} ad_2={ad_2} ad={ad} wd={wd}")

		at = ab + ad*wd
		if self.debug :
			print(f"at={at}")

		i = 4 if ( 0 <= qd * wr ) else 0
		if am < abs(at) :
			self.res[i+1] = [wd, abs(am*wd - ab)/jm]
			self.res[i+2] = [0.0, (at**2 - am**2)/(am*jm)]
			self.res[i+3] = [-wd, abs(am*wd - ab)/jm]
		else :
			self.res[i+1] = [wd, ad / jm]
			self.res[i+3] = [-wd, ad / jm]

		self.sol = list()
		for cmd, dur in self.res :
			print(f"{cmd}\t{dur}")
			if dur > 0.0 :
				self.sol.append([cmd, dur])

		self._integrate_num()

		return self