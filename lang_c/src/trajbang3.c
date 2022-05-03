#ifndef NDEBUG
	#include <stdlib.h>
	#include <stdio.h>
#endif

#include <assert.h>
#include <string.h>

#include "math_lib.h"

#include "trajbang3_mod.h"

static int tb3_discrete_(tb3_C * self, double period, double * jrk);

static int tb3_compute_(tb3_C * self);
static int tb3_integrate_(tb3_C * self);

int tb3_get_qtw(tb3_C * self, double a_from, double a_to, double * q, double * t, double * w) {
	/* return 3 essential parameters for the trajbang computation, to change acceleration, from a_from to a_to at Jm
		q is the area below the curve
		t is the time it takes 
		w is the way the acceleration change
	*/

	double d = a_to - a_from;
	double m = fabs(d);

	*w = (d == 0.0) ? (0.0) : (copysign(1.0, d));

	*t = m / self->jm;
	*q = (a_to + a_from) * (*t) / 2.0;

	return EXIT_SUCCESS;
}

int tb3__init__(tb3_C * self, double jm, double am) {
	// printf("tb3__init__(<self>, jm=%f, am=%f)\n", jm, am);

	memset(self, 0, sizeof(tb3_C));

	self->jm = jm;
	self->am = am;

	return EXIT_SUCCESS;
}

int tb3__setup__(tb3_C * self, double a0, double s0, double ag, double sg) {
	// printf("tb3__setup__(<self>, a0=%f, s0=%f, ag=%f, sg=%f)\n", a0, s0, ag, sg);

	self->a0 = a0;
	self->s0 = s0;
	self->ag = m_BOUND(ag, -self->am, self->am);
	self->sg = sg;

	memset(& (self->cmd), 0, 8 * sizeof(double));
	memset(& (self->dur), 0, 8 * sizeof(double));

	return EXIT_SUCCESS;
}

int tb3__new__(tb3_C * self, double jm, double am, double a0, double s0, double ag, double sg) {
	int are_limits_modified = (jm != self->jm || am != self->am);

	if (are_limits_modified) {
		tb3__init__(self, jm, am);
	}

	int are_conditions_modified = (a0 != self->a0 || s0 != self->s0 || ag != self->ag || sg != self->sg);

	if (are_limits_modified || are_conditions_modified) {
		tb3__setup__(self, a0, s0, ag, sg);

		tb3_compute_(self);
		tb3_integrate_(self);

		assert( tb3_check(self) );
	}

	return EXIT_SUCCESS;
}

static int tb3_compute_(tb3_C * self) {

	double ap = m_BOUND(self->a0, -self->am, self->am);

	double qi, ti, wi;
	tb3_get_qtw(self, self->a0, ap, & qi, & ti, & wi);

	self->cmd[0] = wi;
	self->dur[0] = ti;

	double qr, tr, wr;
	tb3_get_qtw(self, ap, self->ag, & qr, & tr, & wr);

	self->cmd[4] = wr;
	self->dur[4] = tr;

	double qd = self->sg - self->s0 - qi - qr;
	double wd = copysign(1.0, qd);

	double ab = (0 <= qd * wr) ? (self->ag) : (ap);

	double de = ab*ab + self->jm*fabs(qd);
	double ad_1 = ( -ab + sqrt(de) ) * wd;
	double ad_2 = ( -ab - sqrt(de) ) * wd;
	double ad = (0 <= ad_1) ? (ad_1) : (ad_2);	

	double at = ab + ad*wd;

	size_t i = (0 <= qd * wr) ? (4) : (0);
	if ( self->am < fabs(at) ) {
		self->cmd[i+1] = wd;
		self->dur[i+1] = fabs(self->am*wd - ab) / self->jm;

		self->cmd[i+2] = 0.0;
		self->dur[i+2] = (at*at - self->am*self->am) / (self->am*self->jm);

		self->cmd[i+3] = -wd;
		self->dur[i+3] = fabs(self->am*wd - ab) / self->jm;
	} else {
		self->cmd[i+1] = wd;
		self->dur[i+1] = ad / self->jm;

		self->cmd[i+2] = 0.0;
		self->dur[i+2] = 0.0;

		self->cmd[i+3] = -wd;
		self->dur[i+3] = ad / self->jm;
	}

	if (fabs(self->a0 - 2.0) < 1e-6 && fabs(self->s0 - 4.0) < 1e-6) {
		printf("qi=%f, ti=%f, wi=%f\n", qi, ti, wi);
		printf("qr=%f, tr=%f, wr=%f\n", qr, tr, wr);
		printf("ap=%f, ab=%f, at=%f\n", ap, ab, at);
		printf("de=%f, ad_1=%f, ad_1=%f, ad=%f\n", de, ad_1, ad_2, ad);
		for (size_t i=0 ; i<8 ; i++) {
			printf("%f\t%f\n", self->cmd[i], self->dur[i]);
		}
	}

	return EXIT_SUCCESS;
}

int tb3_integrate_(tb3_C * self) {

	for (size_t i=0 ; i < 8 ; i++) {
		self->poly.j[i] = self->cmd[i] * self->jm;
	}

	self->poly.a[0][1] = self->a0;
	self->poly.s[0][2] = self->s0;
	self->poly.p[0][2] = 0.0;

	double t_total = 0.0;
	for (size_t i=0 ; i < 8 ; i++) {
		self->poly.a[i][0] = self->poly.j[i];

		self->poly.s[i][0] = self->poly.a[i][0] / 2.0;
		self->poly.s[i][1] = self->poly.a[i][1];

		self->poly.p[i][0] = self->poly.s[i][0] / 3.0;
		self->poly.p[i][1] = self->poly.s[i][1] / 2.0;
		self->poly.p[i][2] = self->poly.s[i][2];

		double t = self->dur[i];
		t_total += t;
		self->poly.t[i] = t_total;

		self->poly.a[i+1][1] = (
			self->poly.a[i][0]*t +
			self->poly.a[i][1]
		);
		self->poly.s[i+1][2] = (
			self->poly.s[i][0]*t*t + 
			self->poly.s[i][1]*t +
			self->poly.s[i][2]
		);
		self->poly.p[i+1][3] = (
			self->poly.p[i][0]*t*t*t + 
			self->poly.p[i][1]*t*t +
			self->poly.p[i][2]*t +
			self->poly.p[i][3]
		);
	}

	return EXIT_SUCCESS;
}

int tb3_check(tb3_C * self) {

	int error = 0;

	error += m_IS_CLOSE(self->poly.s[8][2], self->sg) ? 0x1 : 0;
	error += m_IS_CLOSE(self->poly.a[8][1], self->ag) ? 0x2 : 0;
	
	for (int i=1 ; i < 9 ; i++) {
		error += m_IS_INSIDE(self->poly.a[i][1], -self->am, true, self->am, true) ? 0x4 << i : 0 ;
	}

	return error;
}

int tb3_at_time(tb3_C * self, double t_req, double * s, double * a, double * j) {

	double t_first = 0.0;
	double t_last = self->poly.t[7];

	if ( m_IS_INSIDE(t_req, t_first, true, t_last, false) ) {
		double t_prev = t_first;
		for ( int i=0 ; i < 8 ; i++ ) {
			double t_next = self->poly.t[i];
			if ( m_IS_INSIDE(t_req, t_prev, true, t_next, false) ) {
				double t = t_req - t_prev;
				*j = self->poly.j[i];
				*a = self->poly.a[i][0]*t + self->poly.a[i][1];
				*s = self->poly.s[i][0]*t*t + self->poly.s[i][1]*t + self->poly.s[i][2];
				return EXIT_SUCCESS;
			}
			t_prev = t_next;
		}
		return 2; // what ? we should have found a solution here
	} else {
		if ( t_req < t_first ) {
			*j = 0.0;
			*a = self->poly.a[0][1];
			*s = self->poly.s[0][2];
			return -1; // time is before the time slot
		}
		if ( t_last <= t_req ) {
			*j = 0.0;
			*a = self->poly.a[8][1];
			*s = self->poly.s[8][2];
			return 1; // time is after the time slot
		}
		return 3; // what ?
	}
}

// int _tb3_compute_1(tb3_C * self, double period); {
// 	// exercice inverse, on accepte de ne pas respecter les limitations, mais en revanche, on cherche à tout faire en 1 step, sans respecter les limites imposées

// 	double qh = period * (self->a0 + self->ag) / 2.0;
// 	double jh = (self->ag - self->a0) / period;

// 	if ( ((self->sg - self.s0) == qh) && ((- self->jm) <= jh) && (jh <= self->jm) ) {
// 		self->cmd[0] = jh;
// 		self->dur[0] = period;

// 		return EXIT_SUCCESS;
// 	}

// 	return EXIT_FAILURE;

// }

static int tb3_discrete_(tb3_C * self, double period, double * jrk) {
	// on cherche à tout faire en 1 ou 2 step, on vérifie les limites après

	double j0 = (( -(self->ag + (3.0 * self->a0)) * period ) + 2.0 * (self->sg - self->s0)) / (2.0 * period * period);
	double j1 = (( ((3.0 * self->ag) + self->a0) * period ) - 2.0 * (self->sg - self->s0)) / (2.0 * period * period);

	double a1 = self->a0 + j0 * period;

	if (m_IS_INSIDE_CLOSED(j0, -self->jm, self->jm) && m_IS_INSIDE_CLOSED(j1, -self->jm, self->jm) && m_IS_INSIDE_CLOSED(a1, -self->am, self->am)) {
		// self->cmd[0] = j0;
		// self->dur[0] = period;

		// self->cmd[1] = j1;
		// self->dur[1] = period;

		* jrk = j0;

		return EXIT_SUCCESS;
	}

	return EXIT_FAILURE;
}

int tb3_jrk(tb3_C * self, double period, double * jrk) {
	
	// tb3__print__(self);

	if (tb3_discrete_(self, period, jrk) == EXIT_SUCCESS) {
		// printf("tb3_discrete() -> %f\n", * jrk);
		return EXIT_SUCCESS;
	}

	double t = 0.0;
	double c = 0.0;

	for (size_t i=0 ; i<8 ; i++) {
		real cmd = self->cmd[i];
		real dur = self->dur[i];
		if ( (t + dur) <= period ) {
			c += cmd * dur;
		} else {
			c += cmd * m_BOUND((period - t), 0.0, period);
		}
		t += dur;
	}

	* jrk = c * self->jm / period;
	// printf("tb3_continuous() -> %f\n", * jrk);
	return EXIT_SUCCESS;

}


int tb3_spd_at_pos(tb3_C * self, double p, double * s) {
	// fprintf(stderr, "----\n");

	if (p <= self->poly.p[0][3]) {
		// fprintf(stderr, "-- ::\n");
		*s = self->poly.s[0][2];
		return EXIT_SUCCESS;
	}
	if (self->poly.p[8][3] <= p) {
		// fprintf(stderr, "++ ::\n");
		*s = self->poly.s[8][2];
		return EXIT_SUCCESS;
	}

	double next_p = self->poly.p[8][3];
	for (int i=7 ; 0 <= i ; i--) {
		double prev_p = self->poly.p[i][3];
		if (p == prev_p) {
			*s = self->poly.s[i][2];
			return EXIT_SUCCESS;
		}
		// fprintf(stderr, "%d :: %.5g < %.5g < %.5g\n", i, prev_p, p, next_p);

		if (prev_p < p && p < next_p) {
			double a = self->poly.p[i][0];
			double b = self->poly.p[i][1];
			double c = self->poly.p[i][2];
			double d = prev_p - p;

			double complex x_lst[3] = {0};
			tb3_poly_sol q_lst[3] = {0};

			for (int j=0 ; j < 3 ; j++) {
				x_lst[j] = NAN;
				q_lst[j] = tb3_poly_sol_NOT_COMPUTED;
			}

			if (a == 0.0) {
				if (b == 0.0) {
					tb3_STATIC_solve_poly1(c, d, x_lst, q_lst);
				} else {
					tb3_STATIC_solve_poly2(b, c, d, x_lst, q_lst);
				}
			} else {
				tb3_STATIC_solve_poly3(a, b, c, d, x_lst, q_lst);
			}

			// for (int j=0 ; j < 3 ; j++) {
			// 	fprintf(stderr, "x%d = %.5g + %.5gi [%d]\n", j, creal(x_lst[j]), cimag(x_lst[j]), q_lst[j]);
			// }

			for (int j=0 ; j < 3 ; j++) {
				if (
					( q_lst[j] == tb3_poly_sol_REAL || q_lst[j] == tb3_poly_sol_MULTIPLE ) &&
					0 <= creal(x_lst[j]) && creal(x_lst[j]) <= self->dur[i]
				) {
					double t = creal(x_lst[j]);
					*s = self->poly.s[i][0]*t*t + self->poly.s[i][1]*t + self->poly.s[i][2];
					return EXIT_SUCCESS;
				}
			}
		}
		next_p = prev_p;
	}

	return EXIT_SUCCESS;
}

int tb3_STATIC_solve_poly1(double a, double b, double complex * x_lst, tb3_poly_sol * q_lst) {

	// solve a x + b = 0

	x_lst[0] = -b / a;
	q_lst[0] = tb3_poly_sol_REAL;

	return EXIT_SUCCESS;
}



int tb3_STATIC_solve_poly2(double a, double b, double c, double complex * x_lst, tb3_poly_sol * q_lst) {

	// solve a x² + b x + c = 0

	// determinant
	double z = b*b - 4.0*a*c;

	if ( 0.0 <= z ) {
		double sqrt_z = sqrt(z);
		x_lst[0] = (- b - sqrt_z)/(2.0*a);
		x_lst[1] = (- b + sqrt_z)/(2.0*a);
		q_lst[0] = q_lst[1] = (z == 0.0) ? (tb3_poly_sol_MULTIPLE) : (tb3_poly_sol_REAL);
	} else {
		double complex csqrt_z = csqrt(z);
		x_lst[0] = (- b - csqrt_z)/(2.0*a);
		x_lst[1] = (- b + csqrt_z)/(2.0*a);
		q_lst[0] = q_lst[1] = tb3_poly_sol_COMPLEX;
	}

	return EXIT_SUCCESS;
}


int tb3_STATIC_solve_poly3_trigo(double a, double b, double c, double d, double complex * x_lst, tb3_poly_sol * q_lst) {
	/* solve a x³ + b x² +cx +d = 0 */
	/* TODO, y a un soucis, un problème de signe qui s'insère pour tb3_STATIC_solve_poly3(a=0.0833 b=0 c=45 d=-55.6) */

	// fprintf(stderr, "tb3_STATIC_solve_poly3(a=%.3g b=%.3g c=%.3g d=%.3g)\n", a, b, c, d);

	// variable change t = x + b / (3 a)
	double e = b/(3.0*a);

	// depressed form t^3 + p t + q = 0
	double p = (3.0*a*c - b*b) / (3.0*a*a);
	double q = (2.0*b*b*b - 9.0*a*b*c + 27.0*a*a*d)/(27.0*a*a*a);

	/* determinants:
		* if z > 0, 3 double solutions
		* if z = 0, 1 double and 2 other double and identical
		* if z < 0, 1 double root and 2 other complex conjugate
	*/
	// double z = - 4.0*p*p*p - 27.0*q*q;

	// for (size_t i=0 ; i < 3 ; i++) {
	// 	x_lst[i] = NAN;
	// 	q_lst[i] = tb3_poly_sol_NOT_COMPUTED;
	// }

	// if ( z <= 0 ) {
	// 	double g = - q / 2.0;
	// 	double h = sqrt( p*p*p/27.0 + q*q/4.0 );
	// 	x_lst[0] = cbrt(g + h) + cbrt(g - h) - e;
	// 	q_lst[0] = tb3_poly_sol_REAL;
	// }

	// fprintf(stderr, "p = %.3g\n", p);
	// fprintf(stderr, "q = %.3g\n", q);
	//// fprintf(stderr, "z = %.3g\n", z);

	double complex s = csqrt(-3.0/p);
	double complex m = ((3.0*q)/(2.0*p))*s;
	// fprintf(stderr, "m = %.3g + %.3g\n", creal(m), cimag(m));
	for (size_t i=0 ; i < 3 ; i++) {
		double complex r = -2.0*s*ccos(cacos(m)/3.0 - 2.0*M_PI*i/3.0)*copysign(1.0, p) - e;
		x_lst[i] = r;
		q_lst[i] = ( m_IS_CLOSE(cimag(r), 0.0) ) ? (tb3_poly_sol_REAL) : (tb3_poly_sol_COMPLEX);
		// // fprintf(stderr, "x%ld = %.5g + %.5gi [%d]\n", i, creal(x_lst[i]), cimag(x_lst[i]), q_lst[i]);
	}

	return EXIT_SUCCESS;
}

int tb3_STATIC_solve_poly3(double a, double b, double c, double d, double complex * x_lst, tb3_poly_sol * q_lst) {
	/* solve a x³ + b x² + cx + d = 0 by the algebraic cardano method
	https://en.wikipedia.org/wiki/Cubic_equation#Cardano's_formula
	*/

	// fprintf(stderr, "tb3_STATIC_solve_poly3(a=%.3g b=%.3g c=%.3g d=%.3g)\n", a, b, c, d);

	// variable change t = x + b / (3 a)
	double e = b/(3.0*a);

	// depressed form t^3 + p t + q = 0
	double p = (3.0*a*c - b*b) / (3.0*a*a);
	double q = (2.0*b*b*b - 9.0*a*b*c + 27.0*a*a*d)/(27.0*a*a*a);

	double complex r = (q*q)/4.0 + (p*p*p)/27.0;
	double complex m = (-q/2.0) + csqrt(r);
	double complex z = cpow(m, 1.0/3.0);

	double complex cbrt_p = (-1 + I*sqrt(3)) * 0.5;
	double complex cbrt_n = (-1 - I*sqrt(3)) * 0.5;

	double complex t_lst[3] = {
		z - p/(3.0*z),
		z*cbrt_p - p/(3.0*z*cbrt_p),
		z*cbrt_n - p/(3.0*z*cbrt_n)
	};

	for (size_t i=0 ; i < 3 ; i++) {
		double complex x = t_lst[i] - e;
		x_lst[i] = x;
		q_lst[i] = ( m_IS_CLOSE(cimag(x), 0.0) ) ? (tb3_poly_sol_REAL) : (tb3_poly_sol_COMPLEX);
		// fprintf(stderr, "x%ld = %.5g + %.5gi [%d]\n", i, creal(x_lst[i]), cimag(x_lst[i]), q_lst[i]);
	}

	return EXIT_SUCCESS;
}

int tb3__print__(tb3_C * self) {
	fprintf(stderr, "jm=%.3g am=%.3g a0=%.3g s0=%.3g ag=%.3g sg=%.3g\n\n", self->jm, self->am, self->a0, self->s0, self->ag, self->sg);

	fprintf(stderr, "dur: ");
	for (size_t i=0 ; i<8 ; i++) {
		fprintf(stderr, "%6.3g |", self->dur[i]);
	}
	fprintf(stderr, "\n");

	fprintf(stderr, "cmd: ");
	for (size_t i=0 ; i<8 ; i++) {
		fprintf(stderr, "%6.3g |", self->cmd[i]);
	}
	fprintf(stderr, "\n\n");

	// fprintf(stderr, "J = \n");
	for (size_t i=0 ; i<8 ; i++) {
		// fprintf(stderr, "\t%ld. %6.3g\n", i, self->poly.j[i]);
	}

	// fprintf(stderr, "A :: af=%.3g ag=%.3g\n", self->poly.a[8][1], self->ag);
	for (size_t i=0 ; i<8 ; i++) {
		// fprintf(stderr, "\t%ld. %6.3g t + %6.3g\n", i, self->poly.a[i][0], self->poly.a[i][1]);
	}

	// fprintf(stderr, "S :: sf=%.3g sg=%.3g\n", self->poly.s[8][2], self->sg);
	for (size_t i=0 ; i<8 ; i++) {
		// fprintf(stderr, "\t%ld. %6.3g t² + %6.3g t + %6.3g\n", i, self->poly.s[i][0], self->poly.s[i][1], self->poly.s[i][2]);
	}

	// fprintf(stderr, "P :: pf=%.3g\n", self->poly.p[8][3]);
	for (size_t i=0 ; i<8 ; i++) {
		// fprintf(stderr, "\t%ld. %6.3g t³ + %6.3g t² + %6.3g t + %6.3g\n", i, self->poly.p[i][0], self->poly.p[i][1], self->poly.p[i][2], self->poly.p[i][3]);
	}

	return EXIT_SUCCESS;
}
