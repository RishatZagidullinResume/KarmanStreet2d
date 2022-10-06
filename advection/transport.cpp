#include "transport.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#define po 1
namespace advection
{

	void inline solve_linear_system_eigen(double *& matrix, double *&result, double *&right_hand_side, int size)
	{
		typedef Eigen::SparseMatrix<double> SpMat;
		typedef Eigen::Triplet<double> T;
		vector<T> tripletList;
		tripletList.reserve(10*size);
		for (int i=0; i<size; i++) for (int j=0; j<size; j++) if (abs(matrix[i*size+j])>1e-4) tripletList.push_back(T(i,j,matrix[i*size+j]));
		SpMat A(size, size);
		Eigen::VectorXd b(size);
		for (int i=0; i<size; i++)
		{
			b[i] = (abs(right_hand_side[i])<1e-4) ? 0.0 : right_hand_side[i];
		}
		A.setFromTriplets(tripletList.begin(), tripletList.end());
		Eigen::LeastSquaresConjugateGradient<SpMat> cg;
		cg.compute(A);
		cg.setTolerance(0.001);
		Eigen::VectorXd x = cg.solve(b);
		//cout << "iterations:	" << cg.iterations() << endl;
		//cout << "estimated error: " << cg.error() << endl;
		for (int i=0; i<size; i++)
		{
			result[i] = x[i];
		}
	}
	

	double advection_2d::limiter(double const &r_factor, double const &weight)
	{
		return max(0.0, min(weight*r_factor,min(0.5*(1.0+r_factor), weight)));//MC
	}

	void advection_2d::find_gradients(Vector2 * values, double * input)
	{
		for (int j = 0; j < size; j++)
		{
			if (!is_boundary[j])
			{	
				double A = 0.0;
				double B = 0.0;
				double C = 0.0;
				double D = 0.0;
				double E = 0.0;
				double F = 0.0;
				for (int k = 0; k < 3; k++)
				{
					if (neighbor_ids[j*3+k] != -1)
					{
						A += pow(barycenter_dist_vector[j*3+k].x(),2)*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
						B += barycenter_dist_vector[j*3+k].x()*barycenter_dist_vector[j*3+k].y()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
						D += pow(barycenter_dist_vector[j*3+k].y(),2)*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);		
						E += barycenter_dist_vector[j*3+k].x()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po)*(input[neighbor_ids[j*3+k]]-input[j]);
						F += barycenter_dist_vector[j*3+k].y()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po)*(input[neighbor_ids[j*3+k]]-input[j]);
					}
				}
				C = B;
				values[j] = Vector2(D/(A*D-B*C)*E - B/(A*D-B*C)*F, - C/(A*D-B*C)*E + A/(A*D-B*C)*F);
			}
			else values[j] = Vector2(0.0, 0.0);
		}
	}

	void advection_2d::calculate_velocities(Vector2 *& velocities, Vector2 *& velocities_original, double *& vorticities, double const & dt, const int& batch)
	{
		for (int i = 0; i < batch; i++)
		{
			double * matrix_coefs = new double [size*size];
			for (int j = 0; j < size*size; j++) matrix_coefs[j] = 0.0;

			double * right_hand_side = new double [size];
			double * result = new double [size];

			for (int j = 0; j < size; j++)
			{

				double A = 0.0;
				double B = 0.0;
				double C = 0.0;
				double D = 0.0;
				double E = 0.0;
				double F = 0.0;
				for (int k = 0; k < 3; k++)
				{
					if (neighbor_ids[j*3+k] != -1)
					{
						A += pow(barycenter_dist_vector[j*3+k].x(),2)*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
						B += barycenter_dist_vector[j*3+k].x()*barycenter_dist_vector[j*3+k].y()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
						D += pow(barycenter_dist_vector[j*3+k].y(),2)*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
						E += -barycenter_dist_vector[j*3+k].x()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
						F += -barycenter_dist_vector[j*3+k].y()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
					}
				}
				C = B;
				//Vector2(D/(A*D-B*C)*E - B/(A*D-B*C)*F, -C/(A*D-B*C)*E + A/(A*D-B*C)*F);
				matrix_coefs[j + size*j] += pow(D/(A*D-B*C)*E -B/(A*D-B*C)*F, 2)  + pow(-C/(A*D-B*C)*E + A/(A*D-B*C)*F, 2);

				for (int k = 0; k < 3; k++)
				{
					if (neighbor_ids[j*3+k] != -1)
					{
						double A_1=0.0;
						double B_1=0.0;
						double C_1=0.0;
						double D_1=0.0;
						double E_1=0.0;
						double F_1=0.0;
						for (int l = 0; l < 3; l++)
						{
							if (neighbor_ids[neighbor_ids[j*3+k]*3+l] != -1)
							{
								A_1 += pow(barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].x(),2)*pow(1.0/barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].length(), po);
								B_1 += barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].x()*barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].y()*pow(1.0/barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].length(), po);
								D_1 += pow(barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].y(),2)*pow(1.0/barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].length(), po);
								E_1 += -barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].x()*pow(1.0/barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].length(), po);
								F_1 += -barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].y()*pow(1.0/barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].length(), po);
							}
						}
						C_1 = B_1;
						double value1 = barycenter_dist_vector[j*3+k].x()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
						double value2 = barycenter_dist_vector[j*3+k].y()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
						matrix_coefs[size*j+neighbor_ids[j*3+k]] += (D/(A*D-B*C)*value1 + -B/(A*D-B*C)*value2)*(D/(A*D-B*C)*E -B/(A*D-B*C)*F)  +  (-C/(A*D-B*C)*value1 + A/(A*D-B*C)*value2)*(-C/(A*D-B*C)*E + A/(A*D-B*C)*F)  +  (D_1/(A_1*D_1-B_1*C_1)*E_1 + -B_1/(A_1*D_1-B_1*C_1)*F_1)*(D/(A*D-B*C)*value1+ -B/(A*D-B*C)*value2) + (-C_1/(A_1*D_1-B_1*C_1)*E_1 + A_1/(A_1*D_1-B_1*C_1)*F_1) * (-C/(A*D-B*C)*value1 + A/(A*D-B*C)*value2);
						for (int l = 0; l < 3; l++)
						{
							if(neighbor_ids[neighbor_ids[j*3+k]*3+l] != -1)
							{
								double another_value1 = barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].x()*pow(1.0/barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].length(), po);
								double another_value2 = barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].y()*pow(1.0/barycenter_dist_vector[neighbor_ids[j*3+k]*3+l].length(), po);
								matrix_coefs[size*j+neighbor_ids[neighbor_ids[j*3+k]*3+l]] += (D/(A*D-B*C)*value1 + -B/(A*D-B*C)*value2)*(D_1/(A_1*D_1-B_1*C_1)*another_value1 + -B_1/(A_1*D_1-B_1*C_1)*another_value2) + (-C/(A*D-B*C)*value1 + A/(A*D-B*C)*value2)*(-C_1/(A_1*D_1-B_1*C_1)*another_value1 + A_1/(A_1*D_1-B_1*C_1)*another_value2);
							}						
						}
					}
				}
				right_hand_side[j] = -vorticities[j];
			}
			solve_linear_system_eigen(matrix_coefs, result, right_hand_side, size);
			find_gradients(&velocities[i*size], result);
			for (int j = 0; j < size; j++)
			{
				if (!is_boundary[j])
				{	
					if (abs(vorticities[j]) >1e-4) velocities[j+size*i] = Vector2(velocities[j+size*i].y(), -velocities[j+size*i].x());
					else velocities[j+size*i] = Vector2(0.0, 0.0); 
				}
				else velocities[j+size*i] = Vector2(0.0, 0.0);
				if (isnan(velocities[j+size*i].x()))
				{
					cout << result[j] << " " << result[neighbor_ids[j*3]] << " "
					     << result[neighbor_ids[j*3+1]] << " " << result[neighbor_ids[j*3+2]] << endl;
					throw runtime_error("isnan(velocities[j*3+k]).x() here");
				}
				if (isnan(velocities[j+size*i].y()))
				{
					throw runtime_error("isnan(face_values[j*3+k]).y()");
				}
			}
			delete [] matrix_coefs;
			delete [] right_hand_side;
			delete [] result;
		}
	}

	void advection_2d::calc_vorticities(Vector2 *& velocities, double *& vorticities, const int& batch)
	{
		for (int i = 0; i < batch; i++)
		{
			double * x_vals = new double [size];
			for (int j = 0; j < size; j++)
			{
				x_vals[j] = velocities[j+size*i].x();
			}
			double * y_vals = new double [size];
			for (int j = 0; j < size; j++)
			{
				y_vals[j] = velocities[j+size*i].y();
			}
			for (int j = 0; j < size; j++)
			{
				int value_index = j+size*i;
				if (!is_boundary[j])
				{	
					double A = 0.0;
					double B = 0.0;
					double C = 0.0;
					double D = 0.0;
					double E = 0.0;
					double F = 0.0;
					double E_2 = 0.0;
					double F_2 = 0.0;
					for (int k = 0; k < 3; k++)
					{
						if (neighbor_ids[j*3+k] != -1)
						{
							A += pow(barycenter_dist_vector[j*3+k].x(),2)*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
							B += barycenter_dist_vector[j*3+k].x()*barycenter_dist_vector[j*3+k].y()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
							D += pow(barycenter_dist_vector[j*3+k].y(),2)*pow(1.0/barycenter_dist_vector[j*3+k].length(), po);
							E += barycenter_dist_vector[j*3+k].x()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po)*(x_vals[neighbor_ids[j*3+k]]-x_vals[j]);
							F += barycenter_dist_vector[j*3+k].y()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po)*(x_vals[neighbor_ids[j*3+k]]-x_vals[j]);
							E_2 += barycenter_dist_vector[j*3+k].x()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po)*(y_vals[neighbor_ids[j*3+k]]-y_vals[j]);
							F_2 += barycenter_dist_vector[j*3+k].y()*pow(1.0/barycenter_dist_vector[j*3+k].length(), po)*(y_vals[neighbor_ids[j*3+k]]-y_vals[j]);
						}
					}
					C = B;
					double val = C/(A*D-B*C)*E - A/(A*D-B*C)*F + D/(A*D-B*C)*E_2 - B/(A*D-B*C)*F_2;
					//if (abs(val-vorticities[value_index])>2.0) vorticities[value_index] += val;
					if (abs(val)>=1e-4) vorticities[value_index] += val;
					//else vorticities[value_index] = 0.0;
					if (isnan(vorticities[value_index]))
					{
						throw runtime_error("isnan(vorticities) first");
					}
				}
				else vorticities[value_index] = 0.0;
			}
			delete [] x_vals;
			delete [] y_vals;
			
		}
	}
	
	void advection_2d::calc_dynamic_arrays(vector<Triangle2*> triangles)
	{
		ofstream barys;
		barys.open("./data/cellcenters.txt");
		is_boundary = new bool [size];
		neighbor_ids = new int [3*size];
		face_normals = new Vector2 [3*size];
		barycenter_dist_vector = new Vector2 [3*size];
		tr_areas = new double [size];
		inverse_weights = new double [3*size];
		is_no_stick = new bool [size];

		int k = 0;
		double leftmost_x = triangles[0]->getCorner(0)->x();
		double rightmost_x = leftmost_x;
		double bottom_y = triangles[0]->getCorner(0)->y();
		double top_y = bottom_y;
		for (auto it(triangles.begin()); it != triangles.end(); it++)
		{
			barys << -(*it)->getBarycenter().y() << " " << (*it)->getBarycenter().x() << " " << 0.0 << endl;
			for (int i = 0; i < 3; i++)
			{
				double value = (*it)->getCorner(i)->x();
				if (value > rightmost_x) rightmost_x = value;
				else if (value < leftmost_x) leftmost_x = value;
				double value_y = (*it)->getCorner(i)->y();
				if (value_y < bottom_y) bottom_y = value_y;
				else if (value_y > top_y) top_y = value_y;
			}
		}
		barys.close();
		for (auto it(triangles.begin()); it != triangles.end(); it++)
		{
			is_boundary[k] = false;
			is_no_stick[k] = false;
			for (int i = 0; i < 3; i++)
			{
				int index = k*3+i;
				auto opposite_tr = (*it)->getOppositeTriangle(i);
				auto neighbor = find(triangles.begin(), triangles.end(), opposite_tr);
				if ((opposite_tr == nullptr) || (neighbor == triangles.end()))
				{
					for (int j=0; j< 3; j++)
					{
						if (((*it)->getCorner(j)->x() == leftmost_x) || ((*it)->getCorner(j)->x() == rightmost_x )/* || (*it)->getCorner(j)->x() <= 0.5*/)
						{
							is_boundary[k] = true;
							break;
						}
						else is_no_stick[k] = true;
					}
					neighbor_ids[index] = -1;
				}
				else
				{
					neighbor_ids[index] = distance(triangles.begin(), neighbor);
				}	
			}
			k++;
		}	
		k = 0;
		for (auto it(triangles.begin()); it != triangles.end(); it++)
		{
			Point2 barry = (*it)->getBarycenter();
			for (int i = 0; i < 3; i++)
			{
				int index = k*3+i;
				
				Point2* corner1 = (*it)->getCorner((i+1)%3);
				Point2* corner2 = (*it)->getCorner((i+2)%3);
				Point2 center = Point2((corner1->x() + corner2->x())/2.0, (corner1->y() + corner2->y())/2.0);
				double length = pow(pow(corner1->x()-corner2->x(), 2) + pow(corner1->y() - corner2->y(), 2), 0.5);
				double x1, y1;
				double x2, y2;
				if (fabs(corner1->y()-corner2->y())<=1e-6)
				{
					x1 = x2 = center.x();
					y1 = center.y() + length;
					y2 = center.y() - length;
				}
				else if (fabs(corner1->x()-corner2->x())<=1e-6)
				{
					y1 = y2 = center.y();
					x1 = center.x() + length;
					x2 = center.x() - length;
				}
				else
				{
					double n_x = (corner1->y()-corner2->y())/(corner2->x()-corner1->x());
					y1 = center.y() + pow(length*length/(1.0 + n_x*n_x),0.5);
					y2 = center.y() - pow(length*length/(1.0 + n_x*n_x),0.5);
					x1 = center.x() + n_x * (y1 - center.y());
					x2 = center.x() + n_x * (y2 - center.y());
				}
				Vector2 dist1 = Vector2(x1 - barry.x(), y1 - barry.y());
				Vector2 dist2 = Vector2(x2 - barry.x(), y2 - barry.y()); 
				Point2 dest = (dist2.length() > dist1.length()) ? Point2(x2, y2) : Point2(x1, y1);
				face_normals[index] = Vector2(dest.x()-center.x(), dest.y()-center.y());
				if (neighbor_ids[index] != -1)
				{
					Point2 neighborbary = (*it)->getOppositeTriangle(i)->getBarycenter();
					inverse_weights[index] = (pow(pow(barry.x()-center.x(),2)+pow(barry.y()-center.y(),2),0.5) + pow(pow(neighborbary.x()-center.x(),2)+pow(neighborbary.y()-center.y(),2),0.5))/pow(pow(barry.x()-center.x(),2)+pow(barry.y()-center.y(),2),0.5);
					//assert(fabs(inverse_weights[index]-2.0) < 1e-3);
				}
				else
				{
					inverse_weights[index] = 0.0;
				}
			}
			k++;
		}
		k = 0;
		for (auto it(triangles.begin()); it != triangles.end(); it++)
		{
			for (int i = 0; i < 3; i++)
			{
				int index = k*3+i;
				if (neighbor_ids[index] != -1)
				{
					Point2 origin = (*it)->getBarycenter();
					Point2 destination = (*it)->getOppositeTriangle(i)->getBarycenter();
					double bary_x, bary_y;
					if (fabs(bary_x = destination.x()-origin.x()) <= 1e-6) bary_x = 0.0;
					if (fabs(bary_y = destination.y()-origin.y()) <= 1e-6) bary_y = 0.0;
					barycenter_dist_vector[index] = Vector2(bary_x, bary_y);
				}
				else barycenter_dist_vector[index] = Vector2(0.0, 0.0);
			}
			k++;
		}
		k = 0;
		for (auto it(triangles.begin()); it != triangles.end(); it++)
		{
			tr_areas[k] = (*it)->getArea2D();
			k++;
		}
	}

	advection_2d::~advection_2d()
	{
		delete [] face_normals;
		delete [] barycenter_dist_vector;
		delete [] neighbor_ids;
		delete [] tr_areas;
		delete [] inverse_weights;
		delete [] is_no_stick;
	}

	void advection_2d::update_velocities(Vector2 *& velocities, const int& batch)
	{
		for (int i = 0; i < batch; i++)
		{
			for (int j = 0; j < size; j++)
			{
				int value_index = j+size*i;
				if (!is_boundary[j])
				{	
					for (int k = 0; k < 3; k++)
					{
						if (neighbor_ids[j*3+k] == -1)
						{
							Vector2 normal = face_normals[j*3+k];
							double coef = velocities[value_index] * normal / (normal.length() * normal.length());
							double x = velocities[value_index].x() - coef * normal.x();
							double y = velocities[value_index].y() - coef * normal.y();
							velocities[value_index] = Vector2(x, y);
							assert(fabs(velocities[value_index]*normal) <= 1e-5);
						}
					}
				}
			}
		}
	}

	void advection_2d::solver(double *&u, double const &dt, Vector2 *&velocities, int const &batch)
	{
		for (int i = 0; i < batch; i++)
		{
			double phi_u_star;

			double * face_values;
			face_values = new double [size * 3];
			Vector2 * cell_gradients;
			cell_gradients = new Vector2 [size];
			double * r_factor_face;
			r_factor_face = new double [size *3];
			Vector2 * interpolated_velocities;
			interpolated_velocities = new Vector2 [size * 3];

			for (int j = 0; j < size; j++)
			{
				int value_index = j+size*i;
				if (!is_boundary[j])
				{	
					for (int k = 0; k < 3; k++)
					{
						if (neighbor_ids[j*3+k] != -1)
						{
							interpolated_velocities[j*3+k] = velocities[neighbor_ids[j*3+k]+size*i]/inverse_weights[j*3+k] + velocities[value_index]/(inverse_weights[j*3+k]/(inverse_weights[j*3+k]-1.0));
						}
						else
						{
							interpolated_velocities[j*3+k] = Vector2(0.0, 0.0);
						}
					}
				}
			}
			find_gradients(cell_gradients, &u[size*i]);
			for (int j = 0; j < size; j++)
			{
				int value_index = j+size*i;
				if (!is_boundary[j])
				{
					for (int k = 0; k < 3; k++)
					{
						double phi_u_star;
						if (neighbor_ids[j*3+k] == -1) face_values[j*3+k] = 0.0;
						else if (face_normals[j*3+k]*interpolated_velocities[j*3+k]<0.0)
						{
							phi_u_star = u[value_index] - 2.0 * (cell_gradients[neighbor_ids[j*3+k]]) * (-barycenter_dist_vector[j*3+k]);
							if (phi_u_star < 0.0) {phi_u_star = 0.0;}
							else if (phi_u_star > 1.0) phi_u_star = 1.0;
							if (fabs(u[value_index] - u[neighbor_ids[j*3+k]+size*i]) > 1e-5)
							{
								r_factor_face[j*3+k] = (u[neighbor_ids[j*3+k]+size*i] - phi_u_star)/(u[value_index]-u[neighbor_ids[j*3+k]+size*i]);
							}
							else
							{
								if (fabs(u[neighbor_ids[j*3+k]+size*i] - phi_u_star) <= 1e-5) r_factor_face[j*3+k] = ((u[neighbor_ids[j*3+k]+size*i] - phi_u_star) >= 0.0) ? 1.0 : -1.0;
								else 
								{
									r_factor_face[j*3+k] = (u[neighbor_ids[j*3+k]+size*i] - phi_u_star) / 1e-5;
								}
							}
							face_values[j*3+k] = u[neighbor_ids[j*3+k]+size*i] + limiter(r_factor_face[j*3+k], inverse_weights[j*3+k]/(inverse_weights[j*3+k]-1.0)) / (inverse_weights[j*3+k]/(inverse_weights[j*3+k]-1.0)) * (u[value_index] - u[neighbor_ids[j*3+k]+size*i]); 
						}
						else
						{
							phi_u_star = u[neighbor_ids[j*3+k]+size*i] - 2.0 * cell_gradients[j] * barycenter_dist_vector[j*3+k];
							if (phi_u_star < 0.0) {phi_u_star = 0.0;}
							else if (phi_u_star > 1.0) phi_u_star = 1.0;
							if (fabs(u[value_index] - u[neighbor_ids[j*3+k]+size*i]) > 1e-5)
							{
								r_factor_face[j*3+k] = (u[value_index] - phi_u_star)/(u[neighbor_ids[j*3+k]+size*i]-u[value_index]);
							}
							else
							{
								if (fabs(u[value_index] - phi_u_star)<= 1e-5) r_factor_face[j*3+k] = ((u[value_index] - phi_u_star) >= 0.0) ? 1.0 : -1.0;
								else r_factor_face[j*3+k] = (u[value_index] - phi_u_star) / 1e-5;
							}
							face_values[j*3+k] = u[value_index] + limiter(r_factor_face[j*3+k], inverse_weights[j*3+k]) / inverse_weights[j*3+k] * (u[neighbor_ids[j*3+k]+size*i] - u[value_index]);
						}
						if (isnan(face_values[j*3+k]))
						{
							throw runtime_error("isnan(face_values[j*3+k]) from solver");
						}
					}
				}
			}
			for (int j = 0; j < size; j++)
			{
				int value_index = j+size*i;
				if (!is_boundary[j])
				{
					double temp = 0.0;
					for (int k = 0; k < 3; k++)
					{
						temp += face_values[k+j*3] * interpolated_velocities[j*3+k] * face_normals[j*3+k]; 
					}
					u[value_index] = u[value_index] - dt / tr_areas[j] * temp;
					//u[value_index] = temp / tr_areas[j];
				}
			}
			delete [] interpolated_velocities;
			delete [] face_values;
			delete [] cell_gradients;
			delete [] r_factor_face;
		}
	}

}
