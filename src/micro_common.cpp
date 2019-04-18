/*
 *  This source code is part of MicroPP: a finite element library
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


#include "micro.hpp"
#include "material.hpp"


template<int tdim>
micropp<tdim>::micropp(const int _ngp, const int size[3], const int _micro_type,
		       const double _micro_params[4],
		       const struct material_base *_materials,
		       const int _coupling, const bool _subiterations,
		       const int _nsubiterations, const int _nr_max_its,
		       const double _nr_max_tol, const double _nr_rel_tol):

	ngp(_ngp),
	nx(size[0]), ny(size[1]),
	nz((tdim == 3) ? size[2] : 1),

	nn(nx * ny * nz),
	nndim(nn * dim),

	nex(nx - 1), ney(ny - 1),
	nez((tdim == 3) ? (nz - 1) : 1),
	coupling(_coupling),

	nelem(nex * ney * nez),
	lx(_micro_params[0]), ly(_micro_params[1]),
	lz((tdim == 3) ? _micro_params[2] : 0.0),
	dx(lx / nex), dy(ly / ney), dz((tdim == 3) ? lz / nez : 0.0),

	special_param(_micro_params[3]),
	subiterations(_subiterations),
	nsubiterations(_nsubiterations),

	wg(((tdim == 3) ? dx * dy * dz : dx * dy) / npe),
	vol_tot((tdim == 3) ? lx * ly * lz : lx * ly),
	ivol(1.0 / (wg * npe)),
	micro_type(_micro_type), nvars(nelem * npe * NUM_VAR_GP),

	nr_max_its(_nr_max_its),
	nr_max_tol(_nr_max_tol),
	nr_rel_tol(_nr_rel_tol)
{
	INST_CONSTRUCT; // Initialize the Intrumentation

	gp_list = new gp_t<tdim>[ngp]();
	for (int gp = 0; gp < ngp; ++gp) {
		gp_list[gp].u_n = (double *) calloc(nndim, sizeof(double));
		gp_list[gp].u_k = (double *) calloc(nndim, sizeof(double));
	}

	elem_type = (int *) calloc(nelem, sizeof(int));
	elem_stress = (double *) calloc(nelem * nvoi, sizeof(double));
	elem_strain = (double *) calloc(nelem * nvoi, sizeof(double));

	int nParams = 4;
	numMaterials = 2;

	for (int i = 0; i < nParams; i++)
		micro_params[i] = _micro_params[i];

	for (int i = 0; i < numMaterials; ++i) {
#ifndef _OPENACC
		material_list[i] = material_t::make_material(_materials[i]);
#else
		material_acc_list[i] = new material_acc(_materials[i]);
#endif
	}

	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {
				const int e_i = glo_elem(ex, ey, ez);
				elem_type[e_i] = get_elem_type(ex, ey, ez);
			}
		}
	}

	memset(ctan_lin, 0.0, nvoi * nvoi * sizeof(double));
	if (coupling != NO_COUPLING)
		calc_ctan_lin();

	for (int gp = 0; gp < ngp; ++gp)
		memcpy(gp_list[gp].ctan, ctan_lin, nvoi * nvoi * sizeof(double));

#ifdef _OPENACC
	ell_init(&A_acc, dim, dim, size);
//#pragma acc enter data copyin(A_acc)
//#pragma acc enter data copyin(A_acc.cols[:A_acc.nrow * A_acc.nnz])
//#pragma acc enter data create(A_acc.vals[:A_acc.nrow * A_acc.nnz], A_acc.r[:A_acc.nrow], A_acc.z[:A_acc.nrow], A_acc.k[:A_acc.nrow], A_acc.p[:A_acc.nrow], A_acc.Ap[:A_acc.nrow])
#endif

}


template <int tdim>
micropp<tdim>::~micropp()
{
	INST_DESTRUCT;

	free(elem_stress);
	free(elem_strain);
	free(elem_type);

	delete [] gp_list;
}


template <int tdim>
int micropp<tdim>::is_non_linear(const int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return (int) gp_list[gp_id].allocated;
}


template <int tdim>
int micropp<tdim>::get_cost(int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return gp_list[gp_id].cost;
}


template <int tdim>
bool micropp<tdim>::has_converged(int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return gp_list[gp_id].converged;
}


template <int tdim>
bool micropp<tdim>::has_subiterated(int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return gp_list[gp_id].subiterated;
}


template <int tdim>
int micropp<tdim>::get_non_linear_gps(void) const
{
	int count = 0;
	for (int gp = 0; gp < ngp; ++gp)
		if (gp_list[gp].allocated)
			count ++;
	return count;
}


template <int tdim>
void micropp<tdim>::calc_ctan_lin()
{

#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < nvoi; ++i) {

		const int ns[3] = { nx, ny, nz };

		ell_matrix A;  // Jacobian
		ell_init(&A, dim, dim, ns);
		double *b = (double *) calloc(nndim, sizeof(double));
		double *du = (double *) calloc(nndim, sizeof(double));
		double *u = (double *) calloc(nndim, sizeof(double));

		double sig_1[6];
		double eps_1[nvoi] = { 0.0 };
		eps_1[i] += D_EPS_CTAN_AVE;

#ifdef _OPENACC
		newton_raphson_acc(&A, b, u, du, eps_1);
#else
		newton_raphson(&A, b, u, du, eps_1);
#endif

		calc_ave_stress(u, sig_1);

		for (int v = 0; v < nvoi; ++v)
			ctan_lin[v * nvoi + i] = sig_1[v] / D_EPS_CTAN_AVE;

		ell_free(&A);
		free(b);
		free(u);
		free(du);
	}
	filter(ctan_lin, nvoi * nvoi, FILTER_REL_TOL);
}


template <int tdim>
material_t *micropp<tdim>::get_material(const int e) const
{
	return material_list[elem_type[e]];
}


template <int tdim>
void micropp<tdim>::get_elem_rhs(const double *u,
				 const double *vars_old,
				 double be[npe * dim],
				 int ex, int ey, int ez) const
{
	INST_START;

	constexpr int npedim = npe * dim;
	double bmat[nvoi][npedim], stress_gp[nvoi], strain_gp[nvoi];

	memset(be, 0, npedim * sizeof(double));

	for (int gp = 0; gp < npe; ++gp) {

		calc_bmat(gp, bmat);

		get_strain(u, gp, strain_gp, ex, ey, ez);
		get_stress(gp, strain_gp, vars_old, stress_gp, ex, ey, ez);

		for (int i = 0; i < npedim; ++i)
			for (int j = 0; j < nvoi; ++j)
				be[i] += bmat[j][i] * stress_gp[j] * wg;
	}
}


template <int tdim>
void micropp<tdim>::get_elem_mat(const double *u,
				 const double *vars_old,
				 double Ae[npe * dim * npe * dim],
				 int ex, int ey, int ez) const
{
	INST_START;
	const int e = glo_elem(ex, ey, ez);
	const material_t *material = get_material(e);

	double ctan[nvoi][nvoi];
	constexpr int npedim = npe * dim;
	constexpr int npedim2 = npedim * npedim;

	double TAe[npedim2] = { 0.0 };

	for (int gp = 0; gp < npe; ++gp) {

		double eps[6];
		get_strain(u, gp, eps, ex, ey, ez);

		const double *vars = (vars_old) ? &vars_old[intvar_ix(e, gp, 0)] : nullptr;
		material->get_ctan(eps, (double *)ctan, vars);

		double bmat[nvoi][npedim], cxb[nvoi][npedim];
		calc_bmat(gp, bmat);

		for (int i = 0; i < nvoi; ++i) {
			for (int j = 0; j < npedim; ++j) {
				double tmp = 0.0;
				for (int k = 0; k < nvoi; ++k)
					tmp += ctan[i][k] * bmat[k][j];
				cxb[i][j] = tmp * wg;
			}
		}

		for (int m = 0; m < nvoi; ++m) {
			for (int i = 0; i < npedim; ++i) {
				const int inpedim = i * npedim;
				const double bmatmi = bmat[m][i];
				for (int j = 0; j < npedim; ++j)
					TAe[inpedim + j] += bmatmi * cxb[m][j];
			}
		}
	}
	memcpy(Ae, TAe, npedim2 * sizeof(double));
}


template <int tdim>
void micropp<tdim>::get_elem_nodes(int n[npe], int ex, int ey, int ez) const
{
	const int nxny = ny * nx;
	const int n0 = ez * nxny + ey * nx + ex;
	n[0] = n0;
	n[1] = n0 + 1;
	n[2] = n0 + nx + 1;
	n[3] = n0 + nx;

	if (dim == 3) {
		n[4] = n[0] + nxny;
		n[5] = n[1] + nxny;
		n[6] = n[2] + nxny;
		n[7] = n[3] + nxny;
	}
}


template<int tdim>
int micropp<tdim>::get_elem_type(int ex, int ey, int ez) const
{
	const double coor[3] = {
		ex * dx + dx / 2.,
		ey * dy + dy / 2.,
		ez * dz + dz / 2. }; // 2D -> dz = 0

	if (micro_type == MIC_SPHERE) { // sphere in the center

		const double rad = special_param;
		const double center[3] = { lx / 2, ly / 2, lz / 2 }; // 2D lz = 0
		double tmp = 0.;
		for (int i = 0; i < dim; ++i)
			tmp += (center[i] - coor[i]) * (center[i] - coor[i]);

		return (tmp < rad * rad);

	} else if (micro_type == MIC_LAYER_Y) { // 2 flat layers in y dir

		const double width = special_param;
		return (coor[1] < width);

	} else if (micro_type == MIC_CILI_FIB_Z) { // a cilindrical fiber in z dir

		const double rad = special_param;
		const double center[3] = { lx / 2, ly / 2, lz / 2 }; // 2D lz = 0
		double tmp = 0.;
		for (int i = 0; i < 2; ++i)
			tmp += (center[i] - coor[i]) * (center[i] - coor[i]);

		return (tmp < rad * rad);

	} else if (micro_type == MIC_CILI_FIB_XZ) { // 2 cilindrical fibers one in x and z dirs

		const double rad = special_param;
		const double cen_1[3] = { lx / 2., ly * .75, lz / 2. };
		double tmp_1 = 0.;
		for (int i = 0; i < 2; ++i)
			tmp_1 += (cen_1[i] - coor[i]) * (cen_1[i] - coor[i]);

		const double cen_2[3] = { lx / 2., ly * .25, lz / 2. };
		double tmp_2 = 0.;
		for (int i = 1; i < 3; ++i)
			tmp_2 += (cen_2[i] - coor[i]) * (cen_2[i] - coor[i]);

		return ((tmp_1 < rad * rad) || (tmp_2 < rad * rad));

	} else if (micro_type == MIC_QUAD_FIB_XYZ) {

	       	/* 3 quad fibers in x, y and z dirs */

		const double width = special_param;
		const double center[3] = { lx / 2., ly / 2., lz / 2. };

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[1] - center[1]) < width)
			return 1;

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[2] - center[2]) < width)
			return 1;

		if (fabs(coor[1] - center[1]) < width &&
		    fabs(coor[2] - center[2]) < width)
			return 1;

		return 0;

	} else if (micro_type == MIC_QUAD_FIB_XZ) {

	       	/* 2 quad fibers in x and z dirs */

		const double width = special_param;
		const double center[3] = { lx / 2., ly / 2., lz / 2. };

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[1] - center[1]) < width)
			return 1;

		if (fabs(coor[1] - center[1]) < width &&
		    fabs(coor[2] - center[2]) < width)
			return 1;

		return 0;

	} else if (micro_type == MIC_QUAD_FIB_XZ_BROKEN_X) {

	       	/* 2 quad fibers in x and z dirs and the one in x is broken */

		const double width = special_param;
		const double center[3] = { lx / 2., ly / 2., lz / 2. };

		if (fabs(coor[0] - center[0]) < width &&
		    fabs(coor[1] - center[1]) < width)
			return 1;

		if (fabs(coor[1] - center[1]) < width &&
		    fabs(coor[2] - center[2]) < width &&
		    (coor[0] < lx * .8 || coor[0] > lx * .9))
			return 1;

		return 0;

	} else if (micro_type == MIC_SPHERES) {

	       	/* Distribution of Several Spheres of diferent sizes */

		const double factor = special_param;

		const int num = 9;
		const double centers[num][3] = {
			{ 0.50 * lx, 0.50 * ly, 0.50 * lz },
			{ 0.22 * lx, 0.65 * ly, 0.10 * lz },
			{ 0.21 * lx, 0.80 * ly, 0.90 * lz },
			{ 0.10 * lx, 0.19 * ly, 0.35 * lz },
			{ 0.20 * lx, 0.17 * ly, 0.85 * lz },
			{ 0.65 * lx, 0.81 * ly, 0.15 * lz },
			{ 0.81 * lx, 0.78 * ly, 0.60 * lz },
			{ 0.77 * lx, 0.35 * ly, 0.25 * lz },
			{ 0.70 * lx, 0.15 * ly, 0.80 * lz }
		};

		double rads[num] = { 0.3, 0.22, 0.13, 0.08, 0.20, 0.17, 0.19, 0.21, 0.15 };
		for (int i = 0; i < num; ++i)
			rads[i] *= factor;

		for (int i = 0; i < num; ++i) {
			double tmp = 0.;
			for (int d = 0; d < dim; ++d)
				tmp += (centers[i][d] - coor[d]) * (centers[i][d] - coor[d]);
			if (tmp < rads[i] * rads[i])
				return 1;
		}

		return 0;
	}

	cerr << "Invalid micro_type = " << micro_type << endl;
	return -1;
}


template <int tdim>
void micropp<tdim>::get_elem_displ(const double *u,
				   double elem_disp[npe * dim],
				   int ex, int ey, int ez) const
{
	int n[npe] ;
	get_elem_nodes(n, ex, ey, ez);

	for (int i = 0 ; i < npe; ++i)
		for (int d = 0; d < dim; ++d)
			elem_disp[i * dim + d] = u[n[i] * dim + d];
}


template <int tdim>
void micropp<tdim>::get_strain(const double *u, int gp, double *strain_gp,
			       int ex, int ey, int ez) const
{
	double elem_disp[npe * dim];
	get_elem_displ(u, elem_disp, ex, ey, ez);

	double bmat[nvoi][npe * dim];
	calc_bmat(gp, bmat);

	memset(strain_gp, 0, nvoi * sizeof(double));
	for (int v = 0; v < nvoi; ++v)
		for (int i = 0; i < npe * dim; i++)
			strain_gp[v] += bmat[v][i] * elem_disp[i];
}


template <int tdim>
void micropp<tdim>::print_info() const
{
	cout << "micropp" << dim << endl;

	cout << "Micro-structure : ";
	switch(micro_type) {
		case(MIC_SPHERE):
			cout << "MIC_SPHERE" << endl;
			break;
		case(MIC_LAYER_Y):
			cout << "MIC_LAYER_Y" << endl;
			break;
		case(MIC_CILI_FIB_Z):
			cout << "MIC_CILI_FIB_Z" << endl;
			break;
		case(MIC_CILI_FIB_XZ):
			cout << "MIC_CILI_FIB_XZ" << endl;
			break;
		case(MIC_QUAD_FIB_XYZ):
			cout << "MIC_QUAD_FIB_XYZ" << endl;
			break;
		case(MIC_QUAD_FIB_XZ):
			cout << "MIC_QUAD_FIB_XZ" << endl;
			break;
		case(MIC_QUAD_FIB_XZ_BROKEN_X):
			cout << "MIC_QUAD_FIB_XZ_BROKEN_X" << endl;
			break;
		case(MIC_SPHERES):
			cout << "MIC_SPHERES" << endl;
			break;
		default:
			cout << "NO TYPE" << endl;
			break;
	}

	cout << "Coupling : ";
	switch(coupling) {
		case(NO_COUPLING):
			cout << "NO_COUPLING" << endl;
			break;
		case(ONE_WAY):
			cout << "ONE_WAY" << endl;
			break;
		case(FULL):
			cout << "FULL" << endl;
			break;
		default:
			break;
	}
       	
	cout << "ngp :" << ngp << " nx :" << nx << " ny :" << ny << " nz :" << nz << " nn :" << nn << endl;
	cout << "lx : " << lx << " ly : " << ly << " lz : " << lz << " param : " << special_param << endl;
	cout << endl;

	for (int i = 0; i < numMaterials; ++i) {
#ifdef _OPENACC
		material_acc_list[i]->print();
#else
		material_list[i]->print();
#endif
		cout << endl;
	}

	cout << "Number of Subiterations :" << nsubiterations << endl;
	cout << endl;
}


template <int tdim>
void micropp<tdim>::get_stress(int gp, const double eps[nvoi],
			       const double *vars_old,
			       double stress_gp[nvoi],
			       int ex, int ey, int ez) const
{
	const int e = glo_elem(ex, ey, ez);
	const material_t *material = get_material(e);
	const double *vars = (vars_old) ? &vars_old[intvar_ix(e, gp, 0)] : nullptr;

	material->get_stress(eps, stress_gp, vars);
}


template <int tdim>
void micropp<tdim>::calc_ave_stress(const double *u, double stress_ave[nvoi],
				    const double *vars_old) const
{
	memset(stress_ave, 0, nvoi * sizeof(double));

	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double stress_aux[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {

					double stress_gp[nvoi], strain_gp[nvoi];
					get_strain(u, gp, strain_gp, ex, ey, ez);
#ifdef _OPENACC
					get_stress_acc(gp, strain_gp, vars_old, stress_gp, ex, ey, ez);
#else
					get_stress(gp, strain_gp, vars_old, stress_gp, ex, ey, ez);
#endif
					for (int v = 0; v < nvoi; ++v)
						stress_aux[v] += stress_gp[v] * wg;

				}
				for (int v = 0; v < nvoi; ++v)
					stress_ave[v] += stress_aux[v];
			}
		}
	}

	for (int v = 0; v < nvoi; ++v)
		stress_ave[v] /= vol_tot;
}


template <int tdim>
void micropp<tdim>::calc_ave_strain(const double *u, double strain_ave[nvoi]) const
{
	memset(strain_ave, 0, nvoi * sizeof(double));

	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double strain_aux[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {
					double strain_gp[nvoi];

					get_strain(u, gp, strain_gp, ex, ey, ez);
					for (int v = 0; v < nvoi; ++v)
						strain_aux[v] += strain_gp[v] * wg;
				}

				for (int v = 0; v < nvoi; v++)
					strain_ave[v] += strain_aux[v];
			}
		}
	}

	for (int v = 0; v < nvoi; v++)
		strain_ave[v] /= vol_tot;
}


template<int tdim>
void micropp<tdim>::calc_fields(double *u, double *vars_old)
{
	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double eps_a[nvoi] = { 0.0 };
				double sig_a[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {

					double stress_gp[nvoi], strain_gp[nvoi];

					get_strain(u, gp, strain_gp, ex, ey, ez);
#ifdef _OPENACC
					get_stress_acc(gp, strain_gp, vars_old, stress_gp, ex, ey, ez);
#else
					get_stress(gp, strain_gp, vars_old, stress_gp, ex, ey, ez);
#endif

					for (int v = 0; v < nvoi; ++v) {
						eps_a[v] += strain_gp[v] * wg;
						sig_a[v] += stress_gp[v] * wg;
					}
				}

				const int e = glo_elem(ex, ey, ez);
				for (int v = 0; v < nvoi; ++v) {
					elem_strain[e * nvoi + v] = eps_a[v] * ivol;
					elem_stress[e * nvoi + v] = sig_a[v] * ivol;
				}
			}
		}
	}
}


/*
 * Evolutes the internal variables for the non-linear material models
 * Calculates the <f_trial_max> max value.
 */

template<int tdim>
bool micropp<tdim>::calc_vars_new(const double *u, const double *_vars_old,
				  double *_vars_new) const
{
	bool non_linear = false;

	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex){

				const int e = glo_elem(ex, ey, ez);
				const material_t *material = get_material(e);

				for (int gp = 0; gp < npe; ++gp) {

					const double *vars_old = (_vars_old) ? &_vars_old[intvar_ix(e, gp, 0)] : nullptr;
					double *vars_new = &_vars_new[intvar_ix(e, gp, 0)];

					double eps[nvoi];
					get_strain(u, gp, eps, ex, ey, ez);

					non_linear |= material->evolute(eps, vars_old, vars_new);
				}
			}
		}
	}

	return non_linear;
}

template<int tdim>
bool micropp<tdim>::calc_vars_new_acc(const double *u, const double *_vars_old,
				      double *_vars_new) const
{
	bool non_linear = false;

	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex){

				const int e = glo_elem(ex, ey, ez);
				const material_acc *material = get_material_acc(e);

				for (int gp = 0; gp < npe; ++gp) {

					const double *vars_old = (_vars_old) ? &_vars_old[intvar_ix(e, gp, 0)] : nullptr;
					double *vars_new = &_vars_new[intvar_ix(e, gp, 0)];

					double eps[nvoi];
					get_strain(u, gp, eps, ex, ey, ez);

					non_linear |= material->evolute(eps, vars_old, vars_new);
				}
			}
		}
	}

	return non_linear;
}


template class micropp<3>;
