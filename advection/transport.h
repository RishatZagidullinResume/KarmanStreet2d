#pragma once

#include <math.h>	
//void solve_linear_system_owl(cusolverSpHandle_t & solver_handle, cusparseHandle_t & handle, cusparseMatDescr_t & descrA, int size, double const * h_A_dense, double const * h_y, double * h_x);

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include "../fade2d/Fade_2D.h"

using namespace GEOM_FADE2D;
using namespace std;

namespace advection
{
	class advection_2d
	{
	protected:
		int size;
		Vector2 * face_normals;
		Vector2 * barycenter_dist_vector;
		double * tr_areas;
		int * neighbor_ids;
		bool * is_no_stick;
		double * inverse_weights;
		double limiter(double const& r_factor, double const& weight);
	public:
		bool * is_boundary;
		//bool * is_no_slip;
	
		void find_gradients(Vector2 * values, double * input);

		void calc_dynamic_arrays(vector<Triangle2*> triangles);
		advection_2d(vector<Triangle2*> triangles) : size(triangles.size()) 
		{
			calc_dynamic_arrays(triangles);
		}
		~advection_2d();
		advection_2d(const advection_2d&) = delete;
		advection_2d& operator= (const advection_2d&) = delete;
		advection_2d& operator= (const advection_2d&&) = delete;
		advection_2d(const advection_2d&&) = delete;
		void update_velocities(Vector2 *& velocities, const int& batch = 1);
		void calc_vorticities(Vector2 *& velocities, double *& vorticities, const int& batch = 1);
		void calculate_velocities(Vector2 *& velocities, Vector2 *& velocities_original, double *& vorticities, double const & dt, const int& batch = 1);

		virtual void solver(double*& values, double const& dt, Vector2*& coefs, const int& batch = 1);
	};
}
