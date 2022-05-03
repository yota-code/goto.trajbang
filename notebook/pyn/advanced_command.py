#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import ipywidgets

import goto.trajbang.saturated_filter as flt

def saturated_filter_demo(epsilon=0.8, omega0=12.0, sat1=3.0, sat2=60.0) :

	step_lst = [0.0,] * 5 + [1.0,] * 995

	ff = flt.Filter2(epsilon, omega0)
	fs = flt.Filter2_sat(epsilon, omega0, sat1, sat2)

	ff_lst = [ff.step(x) for x in step_lst]
	fs_lst = [fs.step(x) for x in step_lst]

	fig = plt.figure()
	axe = fig.add_subplot(1, 1, 1)

	step_line, = plt.plot(step_lst, label="Xcom")
	filtered_line, = plt.plot(ff_lst, label="Xref filtered")
	saturated_line, = plt.plot(fs_lst, label="Xref saturated")

	plt.legend()

	@ipywidgets.interact(
		epsilon=(0.2, 2.0, 0.1), omega0=(1.0, 50.0, 1.0),
		sat1=(0.5, 10.0, 0.1), sat2=(30.0, 100.0, 0.1),
		layout=ipywidgets.widgets.Layout(width='100%')
	)
	def fiddle(epsilon=epsilon, omega0=omega0, sat1=sat1, sat2=sat2) :
		
		ff = flt.Filter2(epsilon, omega0)
		fs = flt.Filter2_sat(epsilon, omega0, sat1, sat2)

		ff_lst = [ff.step(x) for x in step_lst]
		fs_lst = [fs.step(x) for x in step_lst]

		filtered_line.set_ydata(ff_lst)
		saturated_line.set_ydata(fs_lst)

		axe.relim()
		axe.autoscale_view()
		fig.canvas.draw_idle()

	plt.show()


def successive_filter_demo(epsilon=0.8, omega0=12.0, sat1=3.0, sat2=60.0, stop=215) :

	step_lst = [0.0,] * 5 + [1.0,] * 995

	ff = flt.Filter2_sat(epsilon, omega0, sat1, sat2)
	fs = flt.Filter2_sat(0.8, 12.0, 3.0, 60.0)
	fc = flt.Filter2_sat(0.8, 12.0, 3.0, 60.0)

	ff_lst = [ff.step(x) for x in step_lst]
	fc_lst = ff_lst[:stop] + [ff_lst[stop],] * (1000 - stop)
	fs_lst = [fs.step(x) for x in ff_lst]
	fy_lst  = [fc.step(x) for x in fc_lst]

	fig = plt.figure()
	axe = fig.add_subplot(1, 1, 1)

	extracted_line, = plt.plot(fc_lst, '--')
	shy_line, = plt.plot(fy_lst, '--')
	filtered_line, = plt.plot(ff_lst, label="Xcom")
	saturated_line, = plt.plot(fs_lst, label="Xref")

	plt.legend()

	@ipywidgets.interact(
		epsilon=(0.2, 2.0, 0.1), omega0=(1.0, 50.0, 1.0),
		sat1=(0.5, 10.0, 0.1), sat2=(30.0, 100.0, 0.1),
		stop=(10, 990, 1)
	)
	def fiddle(epsilon=epsilon, omega0=omega0, sat1=sat1, sat2=sat2, stop=stop) :
		
		ff = flt.Filter2_sat(epsilon, omega0, sat1, sat2)
		fs = flt.Filter2_sat(0.8, 12.0, 3.0, 60.0)
		fc = flt.Filter2_sat(0.8, 12.0, 3.0, 60.0)

		ff_lst = [ff.step(x) for x in step_lst]
		fc_lst = ff_lst[:stop] + [ff_lst[stop],] * (1000 - stop)
		fs_lst = [fs.step(x) for x in ff_lst]
		fy_lst  = [fc.step(x) for x in fc_lst]

		filtered_line.set_ydata(ff_lst)
		extracted_line.set_ydata(fc_lst)
		saturated_line.set_ydata(fs_lst)
		shy_line.set_ydata(fy_lst)

		axe.relim()
		axe.autoscale_view()
		fig.canvas.draw_idle()

	plt.show()