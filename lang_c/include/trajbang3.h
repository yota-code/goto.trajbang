#ifndef INCLUDE_trajbang3_H
#define INCLUDE_trajbang3_H

typedef struct {

	double j[8];
	double a[8][2];
	double s[8][3];
	double p[8][4];

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

int tb3__init__(tb3_C * self, double jm, double am);
int tb3_get_qtw(tb3_C * self, double a_from, double a_to, double * q, double * t, double * w);
int tb3_compute(tb3_C * self, double a0, double s0, double ag, double sg);
int tb3_integrate(tb3_C * self);
int tb3__print__(tb3_C * self);

#endif /* INCLUDE_trajbang3_H */