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


def symbolic_integration(depth='JSA', length=3) :
	pass

class Tb3Compute() :

	mini_m = 0.001
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
		return self.jm, self.am, self.a0, self.s0, self.ag, self.sg

	def check_is_inside(self, a, low, high) :
		return low - self.check_tolerance <= a <= high + self.check_tolerance

	def check_is_close(self, a, b) :
		return abs(a - b) <= self.check_tolerance

	@property
	def duration(self) :
		return sum(dur for cmd, dur in self.sol)

	@property
	def val(self) :
		return { 'J_m': self.jm, 'A_m': self.am, 'A_0': self.a0, 'S_0': self.s0, 'A_g': self.ag, 'S_g': self.sg, }

	@staticmethod
	def symbolic_expression(length=3) :
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
					a1 * t + a0,
					s2 * t**2 + s1 * t + s0,
					d3 * t**3 + d2 * t**2 + d1 * t + d0
				)

	@property
	def title(self) :
		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg
		return f"TrajBang3(jm={jm}, am={am}, a0={a0}, s0={s0}, ag={ag}, sg={sg})"

	def plot(self, a0, s0, ag, sg, result_dir=None) :
		self.debug = True
		r_comp = self._compute(a0, s0, ag, sg)
		r_inte = self._integrate(a0, s0, ag, sg, r_comp)
		self.debug = False
		print(r_inte)

		t_lst, j_lst, a_lst, s_lst, d_lst = list(), list(), list(), list(), list()
		for i, (c, d) in enumerate(r_comp) :
			t = np.linspace(0, d)
			t_lst.append(t + r_inte['T'][i])

			j_lst.append(np.ones_like(t) * c * self.jm)
			
			a1, a0 = self.poly['A'][i]
			a_lst.append(a1 * t + a0)
			
			s2, s1, s0 = self.poly['S'][i]
			s_lst.append(s2 * t**2 + s1 * t + s0)

			d3, d2, d1, d0 = self.poly['D'][i]
			d_lst.append(d3 * t**3 + d2 * t**2 + d1 * t + d0)
			
		# t_arr = np.hstack(t_lst)
		# j_arr = np.hstack(j_lst)
		# a_arr = np.hstack(a_lst)
		# s_arr = np.hstack(s_lst)
		# d_arr = np.hstack(d_lst)
		
		plt.figure()
		plt.subplot(4, 1, 1)
		for t, j in zip(t_lst, j_lst) :
			plt.step(t, j)
		plt.ylabel('J')
		plt.subplot(4, 1, 2)
		for t, a in zip(t_lst, a_lst) :
			plt.plot(t, a)
		plt.ylabel('A')
		plt.subplot(4, 1, 3)
		for t, s in zip(t_lst, s_lst) :
			plt.plot(t, s)
		plt.ylabel('S')
		plt.subplot(4, 1, 4)
		for t, d in zip(t_lst, d_lst) :
			plt.plot(t, d)
		plt.ylabel('D')
		if result_dir is None :
			plt.show()
		else :
			plt.savefig(result_dir / f"{self.title}.png")		

	def status(self) :
		jm, am, a0, s0, ag, sg = self.jm, self.am, self.a0, self.s0, self.ag, self.sg

		print("cmd :", ' '.join(f"{cmd:7.3g}" for cmd, dur in self.sol))
		print("dur :", ' '.join(f"{dur:7.3g}" for cmd, dur in self.sol))
		# print("-" * 4)
		# print("  T :", ' '.join(f"{i:7.3g}" for i in ([0.0,] + list(itertools.accumulate(T)))))
		# print("  A :", ' '.join(f"{i:7.3g}" for i in A + [ag,]))
		# print("  S :", ' '.join(f"{i:7.3g}" for i in S + [sg,]))
		# print("  D :", ' '.join(f"{i:7.3g}" for i in D))

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

	def run(self, a0, s0, ag, sg) :
		self._prepare(a0, s0, ag, sg)
		self._compute()

	def _compute(self, a0, s0, ag, sg, period=None) :
		if self.debug :
			print(f"Tb3Compute._compute(jm={self.jm}, am={self.am}, a0={a0}, s0={s0}, ag={ag}, sg={sg})")

		r_lst = [[0, 0] for i in range(8)]

		jm, am = self.jm, self.am

		# computation of the initial segment, from A0 to Ap the closest of either Am or -Am
		ap = max(-am, min(a0, am))
		qi, ti, wi = self.get_q(a0, ap)
		r_lst[0] = [wi, ti]

		#computation of the final segment, from Ap to Ag
		qr, tr, wr = self.get_q(ap, ag)
		r_lst[4] = [wr, tr]

		# computation of Qd the remaining part to go from S0 to Sg without the initial and final segments
		qd = sg - s0 - qi - qr
		ab = ag if ( 0 <= qd * wr ) else ap

		# if Qd != 0.0 we need to compute the remaining command
		wd = math.copysign(1.0, qd)

		de = ab**2 + jm*abs(qd)
		ad_1 = ( -ab + math.sqrt(de) ) * wd
		ad_2 = ( -ab - math.sqrt(de) ) * wd
		
		period = None
		period = 0.02

		ad = 0.0
		if (ad_1 < 0.0) and (abs(ad_1) < abs(ad_2)) and (abs(ad_1) < 0.04) :
			ad = ad_1
		if (ad_2 < 0.0) and (abs(ad_2) < abs(ad_1)) and (abs(ad_2) < 0.04) :
			ad = ad_2
		if ad == 0.0 :
			ad = ad_1 if 0.0 <= ad_1 else ad_2

		ad = ad_1 if ab >= 0 else ad_1
			
		at = ab + ad*wd

		i = 4 if ( 0 <= qd * wr ) else 0
		if am < abs(at) :
			print("\n1",  abs(am*wd - ab)/jm, "\n")
			r_lst[i+1] = [wd, abs(am*wd - ab)/jm]
			r_lst[i+2] = [0.0, (at**2 - am**2)/(am*jm)]
			r_lst[i+3] = [-wd, abs(am*wd - ab)/jm]
		else :
			print("\n2",  abs(am*wd - ab)/jm, "\n")
			r_lst[i+1] = [wd, ad / jm]
			r_lst[i+3] = [-wd, ad / jm]

		if self.debug :
			print(f"ap={ap}")
			print(f"qi={qi} ti={ti} wi={wi} (a0 -> ap)")
			print(f"qr={qr} tr={tr} wr={wr} (ap -> ag)")
			print(f"qd={qd}")
			print(f"ab={ab}")
			print(f"ad_1={ad_1} ad_2={ad_2} => ad={ad} wd={wd}")
			print(f"at={at}")
			print([w for w, d in r_lst])
			print([d for w, d in r_lst])
			print("----------------")
		# print("\n===============\n")

		return [[cmd, dur] for cmd, dur in r_lst if 0.0 < dur]

	def _compute_test(self, a0, s0, ag, sg, period=None) :
		print("TESTTESTESTEST")

		if self.debug :
			print(f"Tb3Compute._compute(jm={self.jm}, am={self.am}, a0={a0}, s0={s0}, ag={ag}, sg={sg})")

		r_lst = [[0, 0] for i in range(8)]

		jm, am = self.jm, self.am

		# computation of the initial segment, from A0 to Ap the closest of either Am or -Am
		ap = max(-am, min(a0, am))
		qi, ti, wi = self.get_q(a0, ap)
		r_lst[0] = [wi, ti]

		#computation of the final segment, from Ap to Ag
		qr, tr, wr = self.get_q(ap, ag)
		r_lst[4] = [wr, tr]

		# computation of Qd the remaining part to go from S0 to Sg without the initial and final segments
		qd = sg - s0 - qi - qr
		ab = ag if ( 0 <= qd * wr ) else ap

		# if Qd != 0.0 we need to compute the remaining command
		wd = math.copysign(1.0, qd)

		de = ab**2 + jm*abs(qd)
		ad_1 = ( -ab + math.sqrt(de) ) * wd
		ad_2 = ( -ab - math.sqrt(de) ) * wd
		
		for ad in [ad_1, ad_2] :

			at = ab + ad*wd

			i = 4 if ( 0 <= qd * wr ) else 0
			if am < abs(at) :
				print("\n1",  abs(am*wd - ab)/jm, "\n")
				r_lst[i+1] = [wd, abs(am*wd - ab)/jm]
				r_lst[i+2] = [0.0, (at**2 - am**2)/(am*jm)]
				r_lst[i+3] = [-wd, abs(am*wd - ab)/jm]
			else :
				print("\n2",  abs(am*wd - ab)/jm, "\n")
				r_lst[i+1] = [wd, ad / jm]
				r_lst[i+3] = [-wd, ad / jm]

			if self.debug :
				print(f"ap={ap}")
				print(f"qi={qi} ti={ti} wi={wi} (a0 -> ap)")
				print(f"qr={qr} tr={tr} wr={wr} (ap -> ag)")
				print(f"qd={qd}")
				print(f"ab={ab}")
				print(f"ad_1={ad_1} ad_2={ad_2} => ad={ad} wd={wd}")
				print(f"at={at}")
				print([w for w, d in r_lst])
				print([d for w, d in r_lst])
				print("----------------")
			# print("\n===============\n")

		return [[cmd, dur] for cmd, dur in r_lst if 0.0 < dur]

	def _integrate(self, a0, s0, ag, sg, r_lst) :
		""" numerical integration """

		jm, am = self.jm, self.am
		t0, d0 = 0.0, 0.0

		p_map = { k : list() for k in ['T', 'A', 'S', 'D'] }

		p_map['T'].append(t0)

		print(type(r_lst))
		
		for i, (cmd, dur) in enumerate(r_lst) :
			t = dur
			p_map['T'].append(p_map['T'][-1] + t)

			a1 = cmd * jm
			p_map['A'].append([a1, a0])

			s2 = a1 / 2.0
			s1 = a0
			p_map['S'].append([s2, s1, s0])

			d3 = s2 / 3.0
			d2 = s1 / 2.0
			d1 = s0
			p_map['D'].append([d3, d2, d1, d0])

			a0 = a1 * t + a0
			s0 = s2 * t**2 + s1 * t + s0
			d0 = d3 * t**3 + d2 * t**2 + d1 * t + d0

		p_map['A'].append([a0,])
		p_map['S'].append([s0,])
		p_map['D'].append([d0,])

		self.poly = p_map

		if self.debug :
			for i, (cmd, dur) in enumerate(r_lst) :
				a1, a0 = self.poly['A'][i]
				print(f"A[{i}] = {a1} t + {a0}")

		return p_map

if __name__ == "__main__" :

	u = Tb3Compute(1.0/3.0, 4.0/3.0)
	a0, s0, ag, sg = -0.03333328837510861, 69.99833333558125, -0.06666666666666667, 69.99333333333334
	a0, s0, ag, sg = -0.03333328837510861, 0.99833333558125, -0.06666666666666667, 0.99333333333334
	# a0, s0, ag, sg = 0.0, 0.0, 0.0, 10.0
	#r_lst = u._compute(a0, s0, ag, sg)
	#p_map = u._integrate(a0, s0, ag, sg, r_lst)
	#print(p_map)

	period = 0.1

	j0 = (( -(ag + (3.0 * a0)) * period ) + 2.0 * (sg - s0)) / (2.0 * period**2)
	j1 = (( ((3.0 * ag) + a0) * period ) - 2.0 * (sg - s0)) / (2.0 * period**2)

	a1 = a0 + j0 * period

	print("j0", j0)
	print("j1", j1)
	print("a0", a0)


	# u.plot(a0, s0, ag, sg)