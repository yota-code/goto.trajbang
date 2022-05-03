#!/usr/bin/env python3

"""
draw all possible trajectory for a given set of parameters
"""

import sys

from coordinated_turn import CoordinatedTurn

spd_kt = 100.0

period = 0.01

phi_lim = 20.0
phidot_lim = 16.0
psidot_lim = 8.0
u = CoordinatedTurn(psidot_lim, phi_lim)

psidot_max = u.psidot_deg(spd_kt)
psiddot_max = phidot_lim / phi_lim * psidot_max

print(f"psidot_max = {psidot_max}")
print(f"psiddot_max = {psiddot_max}")

jm = psiddot_max
am = psidot_max

def cmd_iter(am, jm, period) :
	j_lst = list()
	a_lst = [0.0,]
	for i in range(1, 30) :
		j = min( (am - a_lst[-1]) / (jm * period), 1.0 )
		a = a_lst[-1] + j * jm * period
		print(i, (am - sum(j_lst) * jm * period), j_lst, "->", j, sum(j_lst + [j,]) * jm * period)
		j_lst.append(j)
		a_lst.append(a)
		yield a_lst, j_lst

if __name__ =='__main__' :
	import matplotlib.pyplot as plt

	for a_lst, j_lst in cmd_iter(10.0, 1.5, 0.3) :
		u_lst = j_lst

	plt.plot(u_lst)
	plt.show()

	sys.exit(0)

	for j_lst in cmd_iter(am, jm, period) :
		u_lst = j_lst

