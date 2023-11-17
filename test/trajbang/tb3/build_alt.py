#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from cc_pathlib import Path

import goto.trajbang.tb3

speed = 10 # m/s
length = 1000 # m
freq = 10 # Hz

time = length / speed
point_cnt = time * freq


u = goto.trajbang.tb3.Trajectory(0.01, 0.5)
u.compute(0.0, 100.0, 0.0, 200.0)

u.status()

duration = u.duration
print("point_cnt", point_cnt)
print(u.poly['T'])
print("time = ", time)
print("duration = ", u.duration)

z =  time / duration
print(z)
t_arr = np.arange(0.0, u.duration, 1 / (freq * z))

t_lst = list(t_arr) + [u.duration,]
print(len(t_arr), t_lst[-6:])

r_lst = [u.get_at_time(t) for t in t_lst]

a_lst = [r[0] for r in r_lst]
s_lst = [r[1] for r in r_lst]
d_lst = [r[2] for r in r_lst]

plt.figure()
plt.plot([(t * z) for t in t_lst], s_lst)
plt.savefig("degug.png")

header = [
	"#000_ref", "001_ref_dot", "002_mes", "003_mes_dot", "004_value", "005_is_active", "006_mode", "007__I2_alt_fin", "008__I3_distance_to_end"
]

input_lst = [ header, ]
for t, r in zip(t_lst, r_lst) :
	input_lst.append([r[1], r[0]/z, 0.0, 0.0, 0.0, False, 0, 200.0, time - (t*z)])
Path("replay/01/input.tsv").save(input_lst)

input_lst = [ header, ]
for t, r in zip(t_lst, r_lst) :
	if 20.0 <= (t*z) < 40.0 :
		input_lst.append([r[1], r[0]/z, 0.0, 0.0, 160.0, True, 0, 200.0, time - (t*z)])
	else :
		input_lst.append([r[1], r[0]/z, 0.0, 0.0, 160.0, False, 0, 200.0, time - (t*z)])
Path("replay/02/input.tsv").save(input_lst)

input_lst = [ header, ]
for t, r in zip(t_lst, r_lst) :
	if 60.0 <= (t*z) < 80.0 :
		input_lst.append([r[1], r[0]/z, 0.0, 0.0, 140.0, True, 0, 200.0, time - (t*z)])
	else :
		input_lst.append([r[1], r[0]/z, 0.0, 0.0, 140.0, False, 0, 200.0, time - (t*z)])
Path("replay/03/input.tsv").save(input_lst)

input_lst = [ header, ]
for t, r in zip(t_lst, r_lst) :
	if 30.0 <= (t*z) < 60.0 :
		input_lst.append([r[1], r[0]/z, 0.0, 0.0, 120.0, True, 0, 200.0, time - (t*z)])
	else :
		input_lst.append([r[1], r[0]/z, 0.0, 0.0, 140.0, False, 0, 200.0, time - (t*z)])
Path("replay/04/input.tsv").save(input_lst)