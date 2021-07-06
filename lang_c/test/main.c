#include "trajbang3.h"

int main(int argc, char * argv[]) {

	tb3_C u = {0};

	tb3__init__(&u, 1, 2);
	tb3_compute(&u, 0, 10, 0, 0);
	tb3_integrate(&u);

	tb3__print__(&u);

	// double f = 0.0;
	// sscanf(argv[1], "%lf", &f);

	for (int i=0 ; i<100 ; i++) {
		double s = 0;
		double p = (i / 2.0) - 10.0;
		tb3_solve_spd_at_pos(&u, p, &s);
		printf("%f\t%f\n", p, s);
	}
	// double s = 0;
	// tb3_solve_spd_at_pos(&u, f, &s);
	// printf("%f\n", s);

}