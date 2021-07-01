#include "trajbang3.h"

int main(int argc, char * argv[]) {

	tb3_C u = {0};

	tb3__init__(&u, 1, 2);
	tb3_compute(&u, 0, 0, 0, 10);
	tb3_integrate(&u);

	tb3__print__(&u);

}