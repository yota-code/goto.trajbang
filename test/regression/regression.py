#!/usr/bin/env python3

import sys

import numpy as np
import matplotlib.pyplot as plt

from cc_pathlib import Path

import structarray

import goto.trajbang.tb3.compute
import goto.trajbang.tb3.filter

# jm = 1.0
# am = 4.0

# u = goto.trajbang.tb3.compute.Tb3Compute(jm, am, True)

# # at 304.08
# a0, s0, ag, sg = 1.354248584214099, 28.6763344878234, 1.357370414762284, 28.68482283004077


# print( m.run(ag, sg) )

# a0, s0, ag, sg = 1.357370414762284, 28.70345067781317, 1.360551908572638, 28.70686179687501
# m = goto.trajbang.tb3.filter.Tb3Filter(jm, am, a0, s0, period=0.02)
# print( m.run(ag, sg) )


# sys.exit(0)

# r_lst = u._compute(a0, s0-s0, ag, sg-s0)
# print(r_lst)
# print("SUM =", sum(w * d * jm for w, d in r_lst))

# sys.exit(0)

# a0, s0, ag, sg = 1.354248584214099, 28.6763344878234, 1.357370414762284, 28.68482283004077
# u._compute(a0, s0, ag, sg)

# a0, s0, ag, sg = 1.357370414762284, 28.70345067781317, 1.360551908572638, 28.70686179687501
# u._compute(a0, s0, ag, sg)

# a0, s0, ag, sg = 1.357370414762284, 0.70345067781317, 1.360551908572638, 0.70686179687501

# a0, s0, ag, sg = 1.35, 0.0, 1.36, 0.01
# u._compute(a0, s0, ag, sg)
# u.plot(a0, s0-s0, ag, sg-s0)
# a0, s0, ag, sg = 1.354248584214099, 28.6763344878234, 1.357370414762284, 28.68482283004077
# a0, s0, ag, sg = 1.357370414762284, 28.70345067781317, 1.360551908572638, 28.70686179687501

# data_map = {
# 	'20220428_135627@658.68': [-2.124464972071259, 442.0826815195574, -2.12530505027981, 442.0401838193339],
# 	'20220428_135627@658.70': [-2.12530505027981, 442.0401838193339, -2.126962039249389, 442.0378929245877]
# }


# u._compute(* data_map['20220428_135627@658.68'])
# u._compute(* data_map['20220428_135627@658.70'])


record_dir = Path('/dk/vertex/upmv_record')
record_dir = Path('.').resolve()

meta_pth = 'context.tsv'
data_pth = 'clement/upmv_a5457284/20220428_135627.reb'
data_pth = 'clement/upmv_a5457284/20220428_153857.reb'
data_pth = 'clement/upmv_a5457284/20220428_161431.reb'
data_pth = '20220428_174036.reb'

z = structarray.StructArray(record_dir / meta_pth, record_dir / data_pth)
start = 3500

ag_lst = np.array(z['_C_3_upmv_core._C_1_D__root__._C_2_D_ver_override._C_6_tb3_filter._I1_ag'])[start:]
sg_lst = np.array(z['_C_3_upmv_core._C_1_D__root__._C_2_D_ver_override._C_6_tb3_filter._I0_sg'])[start:]

jr_lst = np.array(z['_C_3_upmv_core._C_1_D__root__._C_2_D_ver_override._C_6_tb3_filter._O2_j'])[start:]
ar_lst = np.array(z['_C_3_upmv_core._C_1_D__root__._C_2_D_ver_override._C_6_tb3_filter._O1_a'])[start:]
sr_lst = np.array(z['_C_3_upmv_core._C_1_D__root__._C_2_D_ver_override._C_6_tb3_filter._O0_s'])[start:]

cmd_lst = [
	np.array(z[f'_C_3_upmv_core._C_1_D__root__._C_2_D_ver_override._C_6_tb3_filter.tb3_self.cmd[{i}]'])[start:] for i in range(8)
]

dur_lst = [
	np.array(z[f'_C_3_upmv_core._C_1_D__root__._C_2_D_ver_override._C_6_tb3_filter.tb3_self.dur[{i}]'])[start:] for i in range(8)
]

p_lst = [
	np.array(z[f'_C_3_upmv_core._C_1_D__root__._C_2_D_ver_override._C_6_tb3_filter.tb3_self.{i}'])[start:] for i in ['jm', 'am', 'a0', 's0', 'ag', 'sg']
]

for n, j in enumerate(jr_lst) :
	if j == -1.0 :
		print(n)
		jm, am, a0, s0, ag, sg = [p[n] for p in p_lst]
		u = goto.trajbang.tb3.compute.Tb3Compute(jm, am, True)
		print(">>>>>>>>\n\n\n\n")
		u._compute_test(a0, s0, ag, sg, 0.02)
		print("\n\n\n\n<<<<<<<<<<")
		u._compute(a0, s0, ag, sg, 0.02)
		print([cmd_lst[i][n] for i in range(8)])
		print([dur_lst[i][n] for i in range(8)])
		break

# plt.subplot(2,1,1)
# plt.plot(sg_lst)
# plt.plot(sr_lst)
# plt.subplot(2,1,2)
# plt.plot(ag_lst)
# plt.plot(ar_lst)
# plt.show()

# print(ag_lst - ar_lst)
# sys.exit(0)

#sg_lst[1000:] = 450.0
#ag_lst[1000:] = 0.0

a0 = np.array(z['_C_3_upmv_core._C_1_D__root__._C_2_D_ver_override._C_6_tb3_filter.a_prev'])[start:][0]
s0 = np.array(z['_C_3_upmv_core._C_1_D__root__._C_2_D_ver_override._C_6_tb3_filter.s_prev'])[start:][0]

m = goto.trajbang.tb3.filter.Tb3Filter(jm, am, ag_lst[0], sg_lst[0], period=0.02)

j_lst, af_lst, sf_lst = list(), list(), list()
for a, s, j in zip(ag_lst, sg_lst, jr_lst) :
	# if j == -1.0 :
	# 	jt, at, st = m.run(a, s)
	# 	print("RAAAH", jt, at, st)
	# 	m.obj.debug = True
	# 	m.obj._compute_test(m.a0, m.s0, m.ag, m.sg, m.period)
	# 	m.obj._compute(m.a0, m.s0, m.ag, m.sg, m.period)
	# 	print(m._jerk())
	# 	raise ValueError
	# else :
	
	jt, at, st = m.run(a, s)
	j_lst.append(jt)
	af_lst.append(at)
	sf_lst.append(st)


plt.subplot(3,1,1)
plt.plot(sg_lst)
plt.plot(sr_lst)
plt.plot(sf_lst)
plt.subplot(3,1,2)
plt.plot(ag_lst)
plt.plot(ar_lst)
plt.plot(af_lst)
plt.subplot(3,1,3)
plt.plot(j_lst)
plt.show()