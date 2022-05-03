#ifndef INCLUDE_fctext_trajbang3_mod_H
#define INCLUDE_fctext_trajbang3_mod_H

#include <complex.h>

typedef enum {
	tb3_poly_sol_NOT_COMPUTED,
	tb3_poly_sol_REAL,
	tb3_poly_sol_MULTIPLE,
	tb3_poly_sol_COMPLEX
} tb3_poly_sol;

typedef struct {
	double t[8];

	double j[8];
	
	double a[9][2];
	double s[9][3];
	double p[9][4];
} tb3_poly_T;

typedef struct {
	// maximal values
	double jm;
	double am;
	// initial values
	double a0;
	double s0;
	// target values
	double ag;
	double sg;

	double cmd[8];
	double dur[8];

	tb3_poly_T poly;
} tb3_C;


int tb3__new__(tb3_C * self, double jm, double am, double a0, double s0, double ag, double sg);

int tb3__init__(tb3_C * self, double jm, double am);
int tb3__setup__(tb3_C * self, double a0, double s0, double ag, double sg);

int tb3_get_qtw(tb3_C * self, double a_from, double a_to, double * q, double * t, double * w);

int tb3_check(tb3_C * self);

int tb3_jrk(tb3_C * self, double period, double * jrk);

int tb3_spd_at_pos(tb3_C * self, double p, double * s);
int tb3_at_time(tb3_C * self, double t_req, double * s, double * a, double * j);
int tb3_time_at_pos(tb3_C * self, double p, double * t);

int tb3_STATIC_solve_poly1(double a, double b, double complex * x_lst, tb3_poly_sol * q_lst);
int tb3_STATIC_solve_poly2(double a, double b, double c, double complex * x_lst, tb3_poly_sol * q_lst);
int tb3_STATIC_solve_poly3(double a, double b, double c, double d, double complex * x_lst, tb3_poly_sol * q_lst);

int tb3__print__(tb3_C * self);

#define m_tb3_distance(tb3_obj, distance) { distance = ( (tb3_C *)(& (tb3_obj)) )->poly.p[8][3]; }
#define m_tb3_duration(tb3_obj, duration) { duration = ( (tb3_C *)(& (tb3_obj)) )->poly.t[7]; }

#define m_tb3_spd_at_pos(pos, tb3_obj, spd) ( tb3_spd_at_pos((tb3_C *)(& (tb3_obj)), pos, & spd) )
#define m_tb3_at_time(t, tb3_obj, spd, acc, jrk) ( tb3_at_time((tb3_C *)(& (tb3_obj)), t, & spd, & acc, & jrk) )

#define m_tb3_time_at_pos(p, tb3_self, t) ( tb3_time_at_pos((tb3_C *)(& (tb3_self)), (p), (& (t))) ) 

#define m_tb3_jrk(period, tb3_self, jrk) ( tb3_jrk((tb3_C *)(& tb3_self), period, & jrk) )

#endif /* INCLUDE_fctext_trajbang3_mod_H */