#include "trajbang3.h"

#define MAX(a, b) (a < b ? b : a)
#define MIN(a, b) (a < b ? a : b)
#define IS_CLOSE(a, b, e) (abs(a - b)

int tb3_get_qtw(tb3_C * self, double a_from, double a_to, double * q, double * t, double * w) {

		double d = a_to - a_from;
		double m = fabs(d);

		*w = (d == 0.0) ? 0.0 : copysign(1.0, d);

		*t = m / self->jm;
		*q = (a_to + a_from) * (*t) / 2.0;

		return EXIT_SUCCESS;
}

int tb3__init__(tb3_C * self, double jm, double am) {

	self->jm = jm;
	self->am = am;

	return EXIT_SUCCESS;
}

int tb3_compute(tb3_C * self, double a0, double s0, double ag, double sg) {

	self->a0 = a0;
	self->s0 = s0;
	self->ag = ag;
	self->sg = sg;

	double ap = MAX(-self->am, MIN(self->a0, self->am));

	double qi, ti, wi;
	tb3_get_qtw(self, self->a0, ap, & qi, & ti, & wi);

	self->cmd[0] = wi;
	self->dur[0] = ti;

	double qr, tr, wr;
	tb3_get_qtw(self, ap, self->ag, & qr, & tr, & wr);

	self->cmd[4] = wr;
	self->cmd[4] = tr;

	double qd = self->sg - self->s0 - qi - qr;
	double wd = copysign(1.0, qd);

	double ab = (0 <= qd * qr) ? (self->ag) : (ap);

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

	return EXIT_SUCCESS;

}

int tb3_integrate(tb3_C * self) {

	for (size_t i=0 ; i < 8 ; i++) {
		self->poly.j[i] = self->cmd[i] * self->jm;
	}

	self->poly.a[0][1] = self->a0;
	self->poly.s[0][2] = self->s0;
	self->poly.p[0][2] = 0.0;

	for (size_t i=0 ; i < 8 ; i++) {
		self->poly.a[i][0] = self->poly.j[i];
		//self->poly.a[i][1] = self->af;

		self->poly.s[i][0] = self->poly.a[i][0] / 2.0;
		self->poly.s[i][1] = self->poly.a[i][1];
		//self->poly.s[i][2] = self->sf;

		self->poly.p[i][0] = self->poly.s[i][0] / 3.0;
		self->poly.p[i][1] = self->poly.s[i][1] / 2.0;
		self->poly.p[i][2] = self->poly.s[i][2];
		//self->poly.p[i][3] = self->pf;

		double t = self->dur[i];

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

int tb3_solve_spd_at_pos(tb3_C * self, double p, double * s) {
	fprintf(stderr, "----\n");

	if (p <= self->poly.p[0][3]) {
		fprintf(stderr, "-- ::\n");
		* s = self->poly.s[0][2];
		return EXIT_SUCCESS;
	}
	if (self->poly.p[8][3] <= p) {
		fprintf(stderr, "++ ::\n");
		* s = self->poly.s[8][2];
		return EXIT_SUCCESS;
	}

	double next_p = self->poly.p[8][3];
	for (int i=7 ; 0 <= i ; i--) {
		double prev_p = self->poly.p[i][3];
		if (p == prev_p) {
			* s = self->poly.s[i][2];
			return EXIT_SUCCESS;
		}
		fprintf(stderr, "%d :: %.5g < %.5g < %.5g\n", i, prev_p, p, next_p);

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

			for (int j=0 ; j < 3 ; j++) {
				fprintf(stderr, "x%d = %.5g + %.5gi [%d]\n", j, creal(x_lst[j]), cimag(x_lst[j]), q_lst[j]);
			}

			for (int j=0 ; j < 3 ; j++) {
				if (
					( q_lst[j] == tb3_poly_sol_REAL || q_lst[j] == tb3_poly_sol_MULTIPLE ) &&
					0 <= creal(x_lst[j]) && creal(x_lst[j]) <= self->dur[i]
				) {
					double t = creal(x_lst[j]);
					* s = self->poly.s[i][0]*t*t + self->poly.s[i][1]*t + self->poly.s[i][2];
					return EXIT_SUCCESS;
				}
			}

		}
		next_p = prev_p;
	}
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

int tb3_STATIC_solve_poly3(double a, double b, double c, double d, double complex * x_lst, tb3_poly_sol * q_lst) {
	/* solve a x³ + b x² +cx +d = 0 */

	fprintf(stderr, "tb3_STATIC_solve_poly3(a=%.3g b=%.3g c=%.3g d=%.3g)\n", a, b, c, d);

	// variable change t = x + b / (3 a)
	double e = b/(3.0*a);

	// depressed form t^3 + p t + q = 0
	double p = (3.0*a*c - b*b) / (3.0*a*a);
	double q = (2.0*b*b*b - 9.0*a*b*c + 27.0*a*a*d)/(27.0*a*a*a);

	/* determinants:
		* if z > 0, 3 real solutions
		* if z = 0, 1 real and 2 other real and identical
		* if z < 0, 1 real root and 2 other complex conjugate
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

	fprintf(stderr, "p = %.3g\n", p);
	fprintf(stderr, "q = %.3g\n", q);
	//fprintf(stderr, "z = %.3g\n", z);

	double complex m = ((3.0*q)/(2.0*p))*csqrt(-3.0/p);
	fprintf(stderr, "m = %.3g + %.3g\n", creal(m), cimag(m));
	for (size_t i=0 ; i < 3 ; i++) {
		double complex r = 2.0*csqrt(-p/3.0)*ccos(cacos(m)/3.0 - 2.0*M_PI*i/3.0) - e;
		x_lst[i] = r;
		q_lst[i] = ( cimag(r) == 0.0 ) ? (tb3_poly_sol_REAL) : (tb3_poly_sol_COMPLEX);
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

	fprintf(stderr, "J = \n");
	for (size_t i=0 ; i<8 ; i++) {
		fprintf(stderr, "\t%ld. %6.3g\n", i, self->poly.j[i]);
	}

	fprintf(stderr, "A :: af=%.3g ag=%.3g\n", self->poly.a[8][1], self->ag);
	for (size_t i=0 ; i<8 ; i++) {
		fprintf(stderr, "\t%ld. %6.3g t + %6.3g\n", i, self->poly.a[i][0], self->poly.a[i][1]);
	}

	fprintf(stderr, "S :: sf=%.3g sg=%.3g\n", self->poly.s[8][2], self->sg);
	for (size_t i=0 ; i<8 ; i++) {
		fprintf(stderr, "\t%ld. %6.3g t² + %6.3g t + %6.3g\n", i, self->poly.s[i][0], self->poly.s[i][1], self->poly.s[i][2]);
	}

	fprintf(stderr, "P :: pf=%.3g\n", self->poly.p[8][3]);
	for (size_t i=0 ; i<8 ; i++) {
		fprintf(stderr, "\t%ld. %6.3g t³ + %6.3g t² + %6.3g t + %6.3g\n", i, self->poly.p[i][0], self->poly.p[i][1], self->poly.p[i][2], self->poly.p[i][3]);
	}

	return EXIT_SUCCESS;

}
