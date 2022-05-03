#!/usr/bin/env python3

import sys

import goto.trajbang.tb3.compute
import goto.trajbang.tb3.filter

jm = 1.0
am = 4.0

u = goto.trajbang.tb3.compute.Tb3Compute(jm, am, True)

# # at 304.08
# a0, s0, ag, sg = 1.354248584214099, 28.6763344878234, 1.357370414762284, 28.68482283004077


# m = goto.trajbang.tb3.filter.Tb3Filter(jm, am, a0, s0, period=0.02)
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


a0, s0, ag, sg = 1.354248584214099, 28.6763344878234, 1.357370414762284, 28.68482283004077
u._compute(a0, s0, ag, sg)

a0, s0, ag, sg = 1.357370414762284, 28.70345067781317, 1.360551908572638, 28.70686179687501
u._compute(a0, s0, ag, sg)

