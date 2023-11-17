#!/usr/bin/env python3

import sys

import numpy as np
import matplotlib.pyplot as plt

import sympy

t = sympy.symbols('t')
T1 = sympy.symbols('T1')
Am = sympy.symbols('Am')
Jm = Am / 10.5

T0 = Am / Jm
T2 = T0

j0 = Jm
j1 = 0
j2 = -Jm

a0 = 0
a0_t = j0 * t + a0
a1 = a0_t.subs({t:T0})
a1_t = j1 * t + a1
a2 = a1_t.subs({t:T1})
a2_t = j2 * t + a2
ag = a2_t.subs({t:T2})

print("Ag = ")
sympy.pprint(ag)

s0 = 0
s0_t = sympy.integrate(a0_t, t) + s0
s1 = s0_t.subs({t:T0})
s1_t = sympy.integrate(a1_t, t) + s1
s2 = s1_t.subs({t:T1})
s2_t = sympy.integrate(a2_t, t) + s2
sg = s2_t.subs({t:T2})

print("Sg = ")
sympy.pprint(sg)

sf = sympy.symbols('sf')

T1sol = sympy.solve(sg - sf, T1)[0]
print("--------\nT1 = ")
sympy.pprint(T1sol)

d0 = 0
d0_t = sympy.integrate(s0_t, t) + d0
d1 = d0_t.subs({t:T0})
d1_t = sympy.integrate(s1_t, t) + d1
d2 = d1_t.subs({t:T1})
d2_t = sympy.integrate(s2_t, t) + d2
dg = d2_t.subs({t:T2})

print("Dg = ")
sympy.pprint(dg)
sympy.pprint(dg.subs({T1 : T1sol}))

df = sympy.symbols('df')


Amsol = sympy.solve(dg.subs({T1 : T1sol}) - df, Am)[0]
print("--------\nAmsol = ")
sympy.pprint(Amsol)


sys.exit(0)

T1sol = sympy.solve(sg - sf, T1)[0]

def v(e) :
	try :
		return float(e.subs({T1 : T1sol}).subs({
			Am : float(sys.argv[1]),
			Jm : float(sys.argv[2]),
			sf : float(sys.argv[3]) * 1852 / 3600
		}))
	except :
		print(repr(e))

x_lst = list()
a_lst = list()
s_lst = list()
d_lst = list()

plt.figure()
plt.suptitle(f"time =  {v(T0):.2f} +  {v(T1):.2f} +  {v(T2):.2f} = {v(T0 + T1 + T2):.3f}\nAm = {v(Am)} / Jm = {v(Jm)}")
z = plt.subplot(3,1,1)
plt.plot(np.linspace(0.0, v(T0)), [v(a0_t.subs({t:i})) for i in np.linspace(0.0, v(T0))])
plt.plot(np.linspace(v(T0), v(T1 + T0)), [v(a1_t.subs({t:i})) for i in np.linspace(0, v(T1))])
plt.plot(np.linspace(v(T1 + T0), v(T2 + T1 + T0)), [v(a2_t.subs({t:i})) for i in np.linspace(0, v(T2))])
plt.grid()

plt.subplot(3,1,2, sharex=z)
plt.plot(np.linspace(0.0, v(T0)), [v(s0_t.subs({t:i})) for i in np.linspace(0.0, v(T0))])
plt.plot(np.linspace(v(T0), v(T1 + T0)), [v(s1_t.subs({t:i})) for i in np.linspace(0, v(T1))])
plt.plot(np.linspace(v(T1 + T0), v(T2 + T1 + T0)), [v(s2_t.subs({t:i})) for i in np.linspace(0, v(T2))])
plt.title(f"Sg = {v(sg):.5f}")

plt.grid()

plt.subplot(3,1,3, sharex=z)
plt.plot(np.linspace(0.0, v(T0)), [v(d0_t.subs({t:i})) for i in np.linspace(0.0, v(T0))])
plt.plot(np.linspace(v(T0), v(T1 + T0)), [v(d1_t.subs({t:i})) for i in np.linspace(0, v(T1))])
plt.plot(np.linspace(v(T1 + T0), v(T2 + T1 + T0)), [v(d2_t.subs({t:i})) for i in np.linspace(0, v(T2))])
plt.title(f"Dg = {v(dg):.5f}")
plt.grid()


plt.show()