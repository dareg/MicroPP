/*
 *  This is a test example for MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>

#include <ctime>
#include <cassert>

#include "micro.hpp"

using namespace std;

#define D_EPS 5.0e-4

int main (int argc, char *argv[])
{
	// Execution ./test3d_1 n [print] [steps]
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " <size> " << endl;
		return(1);
	}
	const int n = atoi(argv[1]);
	const int size[3] = { n, n, n };

	const int dir = 0;
	const int time_steps = 10;
	const int micro_type = MIC_CILI_FIB_Z;
	const double micro_params[4] = { 1.0, 1.0, 1.0, 0.2 };

	ofstream file;
	file.open("test_elastic.dat");

	material_base mat_params[2];
	material_set(&mat_params[0], 0, 1.0e6, 0.3, 0.0, 0.0, 0.0);
	material_set(&mat_params[1], 0, 1.0e3, 0.3, 0.0, 0.0, 0.0);

	micropp<3> micro(1, size, micro_type, micro_params, mat_params,
			 NO_COUPLING, true, 5);
	micro.print_info();

	double sig[6];
	double eps[6] = { 0.0 };

	cout << scientific;
	file << scientific << setw(14);
	file << 0.0 << "\t" << 0.0 << endl;

	for (int t = 0; t < time_steps; ++t) {

		cout << "time step = " << t << endl;

		if (t < 80)
			eps[dir] += D_EPS;
		else if (t < 160)
			eps[dir] -= D_EPS;
		else
			eps[dir] += D_EPS;

		micro.set_strain(0, eps);
		micro.homogenize();
		micro.get_stress(0, sig);
		int non_linear = micro.is_non_linear(0);
		int cost = micro.get_cost(0);
		bool has_converged = micro.has_converged(0);

		char filename[128];
		snprintf(filename, 128, "test_elastic_%d", t);
		micro.output(0, filename);

		micro.update_vars();

		cout 	<< "NL        = " << non_linear << endl
			<< "Cost      = " << cost << endl
			<< "Converged = " << has_converged << endl;

		cout << "eps =\t";
		for (int i = 0; i < 6; ++i)
			cout << eps[i] << "\t";
		cout << endl;

		cout << "sig =\t";
		for (int i = 0; i < 6; ++i)
			cout << sig[i] << "\t";
		cout << endl;

		file    << eps[dir] << "\t"
			<< sig[dir] << "\t" << endl;
	}

	file.close();
	return 0;
}
