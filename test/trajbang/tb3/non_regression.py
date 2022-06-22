#!/usr/bin/env python3

import goto.trajbang.tb3

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
	[1, 2, 3, 0, 1, -0.5],
	[1, 4, 1.357370414762284 0.0 1.360551908572638 0.003411119061841106]
]

for jm, am, a0, s0, ag, sg in test_vec :
	u = goto.trajbang.tb3.Trajectory(jm, am, True)
	u.check(a0, s0, ag, sg)