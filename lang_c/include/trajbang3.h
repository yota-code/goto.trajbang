#ifndef INCLUDE_trajbang3_H
#define INCLUDE_trajbang3_H

#include <stdlib.h>
#include <stdio.h>

#include <math.h>

#include <complex.h>

typedef struct {

	double t[9];

	double j[8];

	double a[9][2];
	double s[9][3];
	double p[9][4];

} tb3_poly_T;

typedef struct {

	double jm;
	double am;

	double a0;
	double s0;

	double ag;
	double sg;

	double cmd[8];
	double dur[8];

	tb3_poly_T poly;

} tb3_C;

typedef enum {
	tb3_poly_sol_NOT_COMPUTED,
	tb3_poly_sol_REAL,
	tb3_poly_sol_MULTIPLE,
	tb3_poly_sol_COMPLEX
} tb3_poly_sol;

int tb3__init__(tb3_C * self, double jm, double am);
int tb3_get_qtw(tb3_C * self, double a_from, double a_to, double * q, double * t, double * w);
int tb3_compute(tb3_C * self, double a0, double s0, double ag, double sg);
int tb3_integrate(tb3_C * self);

int tb3_solve_spd_at_pos(tb3_C * self, double p, double * s);

int tb3_STATIC_solve_poly1(double a, double b, double complex * x_lst, tb3_poly_sol * q_lst);
int tb3_STATIC_solve_poly2(double a, double b, double c, double complex * x_lst, tb3_poly_sol * q_lst);
int tb3_STATIC_solve_poly3(double a, double b, double c, double d, double complex * x_lst, tb3_poly_sol * q_lst);

int tb3__print__(tb3_C * self);

#endif /* INCLUDE_trajbang3_H */