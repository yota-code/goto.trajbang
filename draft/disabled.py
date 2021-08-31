#!/usr/bin/env python3

import math
import sys

import numpy as np
import matplotlib.pyplot as plt

import sympy
from sympy import symbols as sym

from cc_pathlib import Path

def sign_var(s) :
	if s == 0.0 :
		return '0'
	else :
		return '+' if 0.0 < s else '-'

def check_target(exp, val) :
	try :
		assert( math.isclose(float(exp), float(val), rel_tol=1e-5, abs_tol=1e-5) )
		return "OK"
	except AssertionError :
		return f"[{exp} # {val}]"

class TrajBang3() :

	def __init__(self, jm, am, a0, s0, sg) :

		# les valeurs max sont toujours positives
		self.jm, self.am, self.a0, self.s0, self.sg = abs(jm), abs(am), a0, s0, sg

		# self.analytical_command()
		# self.analytical_target_pos()

	@property
	def val(self) :
		return {
			'J_m': self.jm,
			'A_m': self.am, 'A_0': self.a0,
			'S_0': self.s0, 'S_g': self.sg,
		}

	def s_float(self, exp) :
		try :
			return float(exp.subs(self.val))
		except AttributeError :
			return float(exp)

	def symbolic_equation(self, n) :
		T = sym(' '.join(f"T_{i}" for i in range(n)))
		J = sym(' '.join(f"J_{i}" for i in range(n)))
		
		A = [sym('A_0'),]
		for i in range(n) :
			A.append( A[-1] + J[i] * T[i] )
			
		S = [sym('S_0'),]
		for i in range(n):
			S.append( S[-1] + A[i] * T[i] + J[i] * T[i]**2 / 2 )
			
		return T, J, A, S

	def target_pos(self) :
		"""
			compute the expected position reached at the end of the trajbang3
		"""
		t = sym('t')

		pos_n = 0
		spd_n = self.s0
		acc_n = self.a0

		for c, d in zip(self.cmd, self.dur) :
			jrk_s = c
			acc_s = sympy.integrate(jrk_s, t) + acc_n
			spd_s = sympy.integrate(acc_s, t) + spd_n
			pos_s = sympy.integrate(spd_s, t) + pos_n
			acc_n = acc_s.subs({'t': d})
			spd_n = spd_s.subs({'t': d})
			pos_n = pos_s.subs({'t': d})
		self.pg = pos_n

		return self.pg

	def compute(self, mode="numeric") :
		if mode == "numeric"  :
			self.cmd, self.dur, self.bch = self._compute_numeric()
		elif mode == "analytic" :
			self.cmd, self.dur, self.bch = self._compute_analytic()
		else :
			raise ValueError

		self.target_pos()

		return self.cmd, self.dur, self.bch		

	def _compute_numeric(self) :
		"""
			This implementation should cover all cases and should be the simplest expression of a valid trajbang3
		"""
		Jm, Am, A0, S0, Sg = self.jm, self.am, self.a0, self.s0, self.sg

		s = math.copysign(1.0, Sg - S0)
		m = 2*Jm*abs(Sg - S0)

		# in 1 step
		Au = math.sqrt(m)
		if A0 == s*Au :
			return [-s], [Au / Jm], "A"

		qm = math.sqrt(2)*Jm

		# in 2 steps
		Tb2 = math.sqrt(A0**2 + m) / qm
		Ta2 = Tb2 - s*A0/Jm
		Av = abs(A0 + s*Ta2*Jm)
		if Av < Am :
			return [s, -s], [Ta2, Tb2], "B"

		z = math.copysign(1.0, A0)

		qn = z*A0**2 / (2*Jm)
		qd = Sg - S0
		qc = qd - qn

		r = math.copysign(1.0, qc * qn)

		# in 3 steps
		Ta = Am - r*abs(A0) / Jm
		Tb = (A0**2 - 2*Am**2 + m) / (2*Am*Jm)
		Tc = (Am / Jm)
		return [z*r, 0, -z*r], [Ta, Tb, Tc], "C"


	def _compute_analytic(self, unfold=False) :
		"""
			this implementation keep symbolic expressions as much as possible, and details all cases
		"""
		bch_lst = list()

		jm, am, a0, s0, sg = self.jm, self.am, self.a0, self.s0, self.sg
		val = self.val

		def u_abs(exp) :
			if unfold :
				p = exp.subs(val)
				if 0 <= p :
					bch_lst.append('+')
					return exp
				else :
					bch_lst.append('-')
					return - exp
			else :
				return abs(exp)

		# creation des symboles
		J_m, A_m, A_0, S_0, S_g = sym('J_m A_m A_0 S_0 S_g')

		# l'intégrale de l'accélération doit être égale au delta de vitesse
		Q_d = S_g - S_0

		if float(A_0.subs(val)) == 0.0 :
			# cas d'un départ à accélération nulle, le sens est donné par le delta de vitesse
			s = 1 if (Q_d).subs(val) >= 0 else -1

			if unfold :
				bch_lst.append(sign_var(s))

			# l'integrale d'une acceleration triangulaire, partant de zero et atteignant l'accéleration maximale
			Q_m = A_m ** 2 / J_m
			if float(abs(Q_d).subs(val)) <= float(Q_m.subs(val)) :
				# BRANCH B
				bch_lst = ["B",] + bch_lst
				# on a pas le temps d'atteindre l'accélération maximale
				# le temps de montée jusqu'à l'accélération requise (pour une aire inférieure ou égale à Q_m)
				T_a = sympy.sqrt(u_abs(Q_d) / J_m)
				cmd = [s, -s]
				dur = [T_a, T_a]
			else :
				# BRANCH A
				bch_lst = ["A",] + bch_lst
				# on a le temps d'atteindre l'accélération maximale, il y a donc une phase plateau à jerk nul
				# time_triangle = am / jm
				T_t = A_m / J_m
				# aire_restante_pour_rectangle = abs(aire_target) - aire_maxi_triangle
				Q_r = u_abs(Q_d) - Q_m
				# time_rectangle = aire_restante_pour_rectangle / am
				T_r = Q_r / A_m
				cmd = [s, 0, -s]
				dur = [T_t, T_r, T_t]
		else :
			s = 1 if a0 >= 0 else -1

			if unfold :
				bch_lst.append(sign_var(s))
				
			# mais pour revenir à zero, il nous faut au moins un triangle, dont l'aire est égale à
			T_n = u_abs(A_0) / J_m
			Q_n = T_n * A_0 / 2
			Q_c = Q_d - Q_n
			if 0 < float((Q_c * Q_n).subs(val)) :
				# l'aire restante est un triangle possiblement tronqué de taille maxi :
				Q_u = u_abs(A_m**2 - A_0**2) / ( J_m )
				if float(abs(Q_c).subs(val)) <= float(Q_u.subs(val)) :
					bch_lst = ["C",] + bch_lst
					D = 4 * A_0**2 + 4 * J_m * u_abs(Q_c)
					T_b = (- 2 * u_abs(A_0) + sympy.sqrt(D) ) / ( 2 * J_m)
					cmd = [s, -s]
					dur = [ T_b, T_b + T_n ]
				else :
					bch_lst = ["D",] + bch_lst
					T_s = (u_abs(Q_c) - Q_u) / A_m
					Q_t = u_abs(A_m - u_abs(A_0))
					cmd = [ s, 0, -s ]
					dur = [ Q_t / J_m, T_s, Q_t / J_m + T_n ]
			elif float((Q_c * Q_n).subs(val)) < 0 :
				# le reste n'est pas dans le même sens que le petit triangle
				# c'est un triangle et pas un trapèze, si l'aire restante est inférieure à :
				Q_s = A_m * A_m / J_m
				if float(Q_s.subs(val)) < float(abs(Q_c).subs(val)) :
					bch_lst = ["E",] + bch_lst
					Q_g = u_abs(Q_c) - Q_s
					cmd = [-s, 0, s]
					dur = [ T_n + A_m / J_m, Q_g / A_m, A_m / J_m ]
				else :
					bch_lst = ["F",] + bch_lst
					T_e = sympy.sqrt(u_abs(Q_c) / J_m)
					cmd = [-s, s ]
					dur = [ T_n + T_e,  T_e]
			else :
				bch_lst = [f"G",] + bch_lst
				cmd = [-s,]
				dur = [T_n,]

		bch = ''.join(bch_lst)

		return cmd, dur, bch

	@property
	def param(self) :
		v = self.val
		f = self.s_float
		return f"Jm={f(v['J_m'])}, Am={f(v['A_m'])}, A0={f(v['A_0'])}, S0={f(v['S_0'])}, Sg={f(v['S_g'])}"

	def prep(self) :
		v = self.val

		def f(exp) :
			try :
				return float(exp.subs(v))
			except AttributeError :
				return float(exp)

		d = [f(i) for i in self.dur] + [1.0,]
		t = np.cumsum([0.0,] + d)
		z = [f(i * v['J_m']) for i in self.cmd] + [0.0,]

		n = len(z)
		for i in range(n) :
			v[f"T_{i}"] = d[i]
			v[f"J_{i}"] = z[i]

		return v, f, n, d, t, z

	def trace(self) :

		v, f, n, d, t, z = self.prep()

		for i in range(n) :
			v[f"T_{i}"] = d[i]
			v[f"J_{i}"] = z[i]
			
		t_xx = [
			np.linspace(t[i], t[i+1], 50)
			for i in range(n)
		]
		j_xx = [
			( np.ones_like(t_xx[0]) * z[i] )
			for i in range(n)
		]
		a_xx = list()
		for i in range(n) :
			a_xx.append(
				(a_xx[-1][-1] if a_xx else f(v['A_0'])) + (t_xx[i] - t[i]) * j_xx[i][0]
			)
		s_xx = list()
		for i in range(n) :
			s_xx.append(
				( s_xx[-1][-1] if s_xx else f(v['S_0']) ) +
				( a_xx[i][0] * (t_xx[i] - t[i]) ) +
				( j_xx[i][0] * (t_xx[i] - t[i])**2 / 2.0 )
			)

		print(f"{self.param} [{self.bch}] cmd={[f(i) for i in self.cmd]} dur={[f(i) for i in self.dur]}")
		if not (
			math.isclose( s_xx[-1][-1], f(v['S_g']), rel_tol=1e-3, abs_tol=1e-3 ) and
			math.isclose( a_xx[-1][-1], 0.0, rel_tol=1e-3, abs_tol=1e-3 )
		) :
			print(f"\terror: Ag={a_xx[-1][-1]:0.3f}, Sg={s_xx[-1][-1]:0.3f}")
			
		self.t_xx, self.j_xx, self.a_xx, self.s_xx = t_xx, j_xx, a_xx, s_xx

	def plot(self, save_dir=None) :
		v, f, n, d, t, z = self.prep()

		t_xx, j_xx, a_xx, s_xx = self.t_xx, self.j_xx, self.a_xx, self.s_xx

		T, J, A, S = self.symbolic_equation(n)

		plt.figure(figsize=(8.27, 11.69))
		title = f"{self.bch} :: {self.param} :: {f(self.pg):0.5e}"

		###  jerk  ###
		plt.subplot(3, 1, 1)
		for i in range(n) :
			plt.plot([t_xx[i][0],] + list(t_xx[i]) + [t_xx[i][-1],], [0.0,] + list(j_xx[i]) + [0.0,])
		plt.grid()
		plt.title(title)
		plt.ylabel("jerk")
		
		###  acceleration  ###
		plt.subplot(3, 1, 2)
		for i in range(n) :
			plt.plot(t_xx[i], a_xx[i])
		print([i.subs(v) for i in A])
		a_lst = [f(i) for i in A]
		plt.plot(t, a_lst, 'o', color="grey")
		plt.grid()
		plt.ylabel("acceleration")
		
		###  speed  ###
		plt.subplot(3, 1, 3)
		for i in range(n) :
			plt.plot(t_xx[i], s_xx[i])
		s_lst = [f(i) for i in S]
		plt.plot(t, s_lst, 'o', color="grey")
		plt.grid()
		plt.ylabel("speed")

		# check_sg, check_ag = check_target(s_lst[-1], val['S_g']), check_target(a_lst[-1], 0.0)

		if save_dir :
			save_pth = Path(save_dir) / f"trajbang3({self.param}).png"
			Path(save_dir).make_dirs()
			plt.savefig(str(save_pth))
			print(save_pth)
			plt.close()
		else :
			save_pth = None
			plt.show()	

