/*
 * This source code is part of MicroPP: a finite element library
 * to solve microstructural problems for composite materials.
 *
 * Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                        Guido Giuntoli <gagiuntoli@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include "material_base.h"

#include <cmath>
#include <cstdio>

struct material_t : public material_base {

	material_t() {
		E = NAN;
		nu = NAN;
		Ka = NAN;
		Sy = NAN;
		k = NAN;
		mu = NAN;
		lambda = NAN;
		type = -1;
		plasticity = false;
		damage = false;
	}

	void set(double _E, double _nu, double _Ka, double _Sy, int _type)
	{
		material_set(this, _E, _nu, _Ka, _Sy, _type);
	}

	void print() const
	{
		material_print(this);
	}
};


#endif
