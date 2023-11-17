#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import signal

NE MARCHE PAS TEL QUEL

dt = 10 # ms

s_arr = np.array([0.0,] * 10 + [1.0] * 80 + [0.5] * 80)

srl = signal.RateLimiter(1.0)
arl = signal.RateLimiter(8.0)
der = signal.Derivative()
ite = signal.Integrate()

pem = signal.JerkLimiter(1.0, 8.0)

z_arr = np.array([srl.run(x) for x in s_arr])
q_arr = np.array([der.run(x) for x in z_arr])

a_arr = np.array([arl.run(x) for x in q_arr])
p_arr = np.array([ite.run(x) for x in a_arr])

m_arr = np.array([pem.run(x) for x in s_arr])

plt.subplot(2,1,1)
plt.plot(s_arr)
plt.plot(z_arr)
plt.plot(p_arr)
plt.plot(m_arr)
plt.grid()
plt.subplot(2,1,2)
plt.plot(q_arr)
plt.plot(a_arr)
plt.grid()

plt.show()