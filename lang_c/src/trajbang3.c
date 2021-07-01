#include <stdlib.h>
#include <stdio.h>

#include <math.h>

#include <trajbang3.h>

#define MAX(a, b) (a < b ? b : a)
#define MIN(a, b) (a < b ? a : b)

int tb3_get_qtw(tb3_C * self, double a_from, double a_to, double * q, double * t, double * w) {

		double d = a_to - a_from;
		double m = abs(d);

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

	double de = ab*ab + self->jm*abs(qd);
	double ad_1 = ( -ab + sqrt(de) ) * wd;
	double ad_2 = ( -ab - sqrt(de) ) * wd;
	double ad = (0 <= ad_1) ? (ad_1) : (ad_2);	

	double at = ab + ad*wd;

	size_t i = (0 <= qd * wr) ? (4) : (0);
	if ( self->am < abs(at) ) {
		self->cmd[i+1] = wd;
		self->dur[i+1] = abs(self->am*wd - ab) / self->jm;

		self->cmd[i+2] = 0.0;
		self->dur[i+2] = (at*at - self->am*self->am) / (self->am*self->jm);

		self->cmd[i+3] = -wd;
		self->dur[i+3] = abs(self->am*wd - ab) / self->jm;
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

	double prev_a = 0.0;
	for (size_t i=0 ; i < 8 ; i++) {
		self->poly.a[i][0] = self->poly.j[i];
		self->poly.a[i][1] = prev_a;
		prev_a = (
			self->poly.a[i][0]*self->dur[i] +
			self->poly.a[i][1]
		);
	}

	double prev_s = 0.0;
	for (size_t i=0 ; i < 8 ; i++) {
		self->poly.s[i][0] = self->poly.a[i][0] / 2.0;
		self->poly.s[i][1] = self->poly.a[i][1];
		self->poly.s[i][2] = prev_s;
		prev_s = (
			self->poly.s[i][0]*self->dur[i]*self->dur[i] + 
			self->poly.s[i][1]*self->dur[i] +
			self->poly.s[i][2]
		);
	}

	double prev_p = 0.0;
	for (size_t i=0 ; i < 8 ; i++) {
		self->poly.p[i][0] = self->poly.s[i][0] / 3.0;
		self->poly.p[i][1] = self->poly.s[i][1] / 2.0;
		self->poly.p[i][2] = self->poly.s[i][2];
		self->poly.p[i][3] = prev_p;
		double t = self->dur[i];
		prev_p = (
			self->poly.p[i][0]*t*t*t + 
			self->poly.p[i][1]*t*t +
			self->poly.p[i][2]*t +
			self->poly.p[i][3]
		);
	}

}

int tb3__print__(tb3_C * self) {

	printf("jm=%f am=%f a0=%f s0=%f ag=%f sg=%f\n", self->jm, self->am, self->a0, self->s0, self->ag, self->sg);

	printf("dur: ");
	for (size_t i=0 ; i<8 ; i++) {
		printf("%6.3g |", self->dur[i]);
	}
	printf("\n");

	printf("cmd: ");
	for (size_t i=0 ; i<8 ; i++) {
		printf("%6.3g |", self->cmd[i]);
	}
	printf("\n\n");

	printf("J = \n");
	for (size_t i=0 ; i<8 ; i++) {
		printf("\t%ld. %6.3g\n", i, self->poly.j[i]);
	}

	printf("A = \n");
	for (size_t i=0 ; i<8 ; i++) {
		printf("\t%ld. %6.3g t + %6.3g\n", i, self->poly.a[i][0], self->poly.a[i][1]);
	}

	printf("S = \n");
	for (size_t i=0 ; i<8 ; i++) {
		printf("\t%ld. %6.3g t² + %6.3g t + %6.3g\n", i, self->poly.s[i][0], self->poly.s[i][1], self->poly.s[i][2]);
	}

	printf("P = \n");
	for (size_t i=0 ; i<8 ; i++) {
		printf("\t%ld. %6.3g t³ + %6.3g t² + %6.3g t + %6.3g\n", i, self->poly.p[i][0], self->poly.p[i][1], self->poly.p[i][2], self->poly.p[i][3]);
	}

	return EXIT_SUCCESS;

}
