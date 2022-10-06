#include <algorithm>
#include <memory>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "include_fade2d/Fade_2D.h"
#include <cassert>
using namespace GEOM_FADE2D;
using namespace std;

class advection_include_2d
{
private:
	double dt;
	int Xrange;
	bool * is_tr_boundary;
	Vector2 * velocities;
	Vector2 * face_normals;
	Vector2 * barycenter_dist_vector;
	Vector2 * skewness_corrector;
	int * neighbor_ids;
	double * tr_areas;
	double * L;
	void calc_is_tr_boundary(vector<Triangle2 *> &triangles);
	void calc_vel(vector<Triangle2 *> &triangles);
	void calc_face_normals_and_more(vector<Triangle2 *> &triangles);
	void calc_bary_dist_vec(vector<Triangle2 *> &triangles);
	void calc_neighbor_ids(vector<Triangle2 *> &triangles);
	void calc_tr_areas(vector<Triangle2 *> &triangles);
	double limiter(double r_factor, double L, double value_downwind, double value_center, Vector2 face_gradient, Vector2 distance);
public:
	double * u;
	advection_include_2d(vector<Triangle2 *> &triangles, double dt);
	~advection_include_2d();
	void solver(int t);
};

advection_include_2d::advection_include_2d(vector<Triangle2 *> &triangles, double dt)
{
	this->dt = dt;
	this->Xrange = triangles.size();
	u = new double [Xrange];
	L = new double [Xrange*3];
	is_tr_boundary = new bool [Xrange];
	velocities = new Vector2 [Xrange];
	face_normals = new Vector2 [Xrange*3];
	skewness_corrector = new Vector2 [Xrange*3];
	barycenter_dist_vector = new Vector2 [Xrange*3];
	neighbor_ids = new int [Xrange*3];
	tr_areas = new double [Xrange];
	calc_is_tr_boundary(triangles);
	calc_neighbor_ids(triangles);
	calc_vel(triangles);	
	calc_face_normals_and_more(triangles);
	calc_bary_dist_vec(triangles);
	calc_tr_areas(triangles);
	
	int k = 0;
	for (auto it(triangles.begin()); it != triangles.end(); it++)
	{
		
		//u[k] = exp(-pow((*it)->getBarycenter().y()-0.5,2)*200.0-pow((*it)->getBarycenter().x()-0.2,2)*200.0);
		//if (((*it)->getBarycenter().x() < 0.3) && ((*it)->getBarycenter().x() > 0.1) && ((*it)->getBarycenter().y() < 0.6) && ((*it)->getBarycenter().y() > 0.4)) u[k] = 1.0;
		//else u[k] = 0.0;
		if (is_tr_boundary[k])
		{
			int how_many = 0;
			Point2 * point;
			for (int i = 0; i < 3; i++)
			{
				point = (*it)->getCorner(i);
				if (point->x() == 0.0) how_many ++;
			}
			if (how_many >= 1) u[k] = 1.0;//exp(-25.0*pow((*it)->getBarycenter().y()-0.5,2));
			else u[k] = 0.0;
		}
		else
		{
			u[k] = 0.0;
		}
		k++;
	}
}

advection_include_2d::~advection_include_2d()
{
	delete [] is_tr_boundary;
	delete [] velocities;
	delete [] face_normals;
	delete [] barycenter_dist_vector;
	delete [] neighbor_ids;
	delete [] tr_areas;
	delete [] skewness_corrector;
	delete [] u;
	delete [] L;
}

void advection_include_2d::solver(int t)
{	
	double * face_values;
	face_values = new double [Xrange * 3];
	Vector2 * cell_gradients;
	cell_gradients = new Vector2 [Xrange];
	double * r_factor_face;
	r_factor_face = new double [Xrange *3];
	for (int j = 0; j < Xrange; j++)
	{
		if (!is_tr_boundary[j])
		{	
			double A, B, C, D, E, F;
			A = pow(barycenter_dist_vector[j*3+0].x()*1.0/barycenter_dist_vector[j*3+0].length(), 2) + pow(barycenter_dist_vector[j*3+1].x()*1.0/barycenter_dist_vector[j*3+1].length(), 2) + pow(barycenter_dist_vector[j*3+2].x()*1.0/barycenter_dist_vector[j*3+2].length(), 2);
			B = C = barycenter_dist_vector[j*3+0].x()*barycenter_dist_vector[j*3+0].y()*pow(1.0/barycenter_dist_vector[j*3+0].length(), 2) + barycenter_dist_vector[j*3+1].x()*barycenter_dist_vector[j*3+1].y()*pow(1.0/barycenter_dist_vector[j*3+1].length(), 2) + barycenter_dist_vector[j*3+2].x()*barycenter_dist_vector[j*3+2].y()*pow(1.0/barycenter_dist_vector[j*3+2].length(), 2);
			D = pow(barycenter_dist_vector[j*3+0].y()*1.0/barycenter_dist_vector[j*3+0].length(), 2) + pow(barycenter_dist_vector[j*3+1].y()*1.0/barycenter_dist_vector[j*3+1].length(), 2) + pow(barycenter_dist_vector[j*3+2].y()*1.0/barycenter_dist_vector[j*3+2].length(), 2);;
			E = barycenter_dist_vector[j*3+0].x()*pow(1.0/barycenter_dist_vector[j*3+0].length(), 2)*(u[neighbor_ids[j*3+0]]-u[j]) + barycenter_dist_vector[j*3+1].x()*pow(1.0/barycenter_dist_vector[j*3+1].length(), 2)*(u[neighbor_ids[j*3+1]]-u[j]) + barycenter_dist_vector[j*3+2].x()*pow(1.0/barycenter_dist_vector[j*3+2].length(), 2)*(u[neighbor_ids[j*3+2]]-u[j]);
			F = barycenter_dist_vector[j*3+0].y()*pow(1.0/barycenter_dist_vector[j*3+0].length(), 2)*(u[neighbor_ids[j*3+0]]-u[j]) + barycenter_dist_vector[j*3+1].y()*pow(1.0/barycenter_dist_vector[j*3+1].length(), 2)*(u[neighbor_ids[j*3+1]]-u[j]) + barycenter_dist_vector[j*3+2].y()*pow(1.0/barycenter_dist_vector[j*3+2].length(), 2)*(u[neighbor_ids[j*3+2]]-u[j]);
			cell_gradients[j] = Vector2(D/(A*D-B*C)*E-B/(A*D-B*C)*F, -C/(A*D-B*C)*E + A/(A*D-B*C)*F);
			//cell_gradients[j] = Vector2(0.0, 0.0);
			//for (int k = 0; k < 3; k++)
			//{
				//face_values[j*3+k] = u[j]*(1.0-1.0/L[j*3+k]) + u[neighbor_ids[j*3+k]]*1.0/L[j*3+k];
				//face_values[j*3+k] = 1.0;
				//face_values[j*3+k] = u[j]*(1.0/2.0) + u[neighbor_ids[j*3+k]]*1.0/2.0;
				
				//cell_gradients[j] = cell_gradients[j] + face_values[j*3+k]*face_normals[j*3+k];
			//}
			//cell_gradients[j] = cell_gradients[j] / tr_areas[j];
		}
		//else cell_gradients[j] = Vector2(0.0, 0.0);
	}
	for (int j = 0; j < Xrange; j++)
	{
		if (!is_tr_boundary[j])
		{
			for (int k = 0; k < 3; k++)
			{
				double phi_u_star;
				//assert((face_normals[j*3+k]*velocities[j])*(barycenter_dist_vector[j*3+k]*velocities[j])>=0.0);
				if (barycenter_dist_vector[j*3+k]*velocities[j]<0.0)
				{
					phi_u_star = u[j] - 2.0 * cell_gradients[neighbor_ids[j*3+k]] * (-barycenter_dist_vector[j*3+k]);
					if (phi_u_star < 0.0) phi_u_star = 0.0;
					if (u[j] != u[neighbor_ids[j*3+k]])
					{
						//r_factor_face[j*3+k] = (2.0 * cell_gradients[neighbor_ids[j*3+k]] * (-barycenter_dist_vector[j*3+k]))/(u[j]-u[neighbor_ids[j*3+k]])-1.0;
						r_factor_face[j*3+k] = (u[neighbor_ids[j*3+k]] - phi_u_star)/(u[j]-u[neighbor_ids[j*3+k]]);
					}
					else
					{
						//cout << "else is here :/\n";
						//r_factor_face[j*3+k] = (2.0 * cell_gradients[neighbor_ids[j*3+k]] * (-barycenter_dist_vector[j*3+k]))/(1e-40)-1.0;
						if (fabs(u[neighbor_ids[j*3+k]] - phi_u_star) <= 1e-6) r_factor_face[j*3+k] = ((u[neighbor_ids[j*3+k]] - phi_u_star) >= 0.0) ? 1.0 : -1.0;
						else r_factor_face[j*3+k] = (u[neighbor_ids[j*3+k]] - phi_u_star)/(1e-11);
						//cout << u[neighbor_ids[j*3+k]] - phi_u_star << endl;
					}
					//face_values[j*3+k] = u[neighbor_ids[j*3+k]] + 0.5 * limiter(r_factor_face[j*3+k], L[j*3+k]) * (u[j]-u[neighbor_ids[j*3+k]]);
					//face_values[j*3+k] = u[neighbor_ids[j*3+k]] + limiter(r_factor_face[j*3+k], L[neighbor_ids[j*3+k]])/L[neighbor_ids[j*3+k]] * (u[j] - u[neighbor_ids[j*3+k]]); 
					face_values[j*3+k] = u[neighbor_ids[j*3+k]] + limiter(r_factor_face[j*3+k], (L[j*3+k]/(L[j*3+k]-1.0)), u[j], u[neighbor_ids[j*3+k]], (u[j]*(1.0-1.0/L[j*3+k]) + u[neighbor_ids[j*3+k]]*1.0/L[j*3+k])*(-face_normals[j*3+k]), skewness_corrector[j*3+k])/(L[j*3+k]/(L[j*3+k]-1.0)) * (u[j] - u[neighbor_ids[j*3+k]]); 
				}
				else// if (barycenter_dist_vector[j*3+k]*velocities[j]>0.0)
				{
					phi_u_star = u[neighbor_ids[j*3+k]] - 2.0 * cell_gradients[j] * barycenter_dist_vector[j*3+k];
					if (phi_u_star < 0.0) phi_u_star = 0.0;
					if (u[j] != u[neighbor_ids[j*3+k]])
					{
						//r_factor_face[j*3+k] = (2.0 * cell_gradients[j] * barycenter_dist_vector[j*3+k])/(u[neighbor_ids[j*3+k]]-u[j])-1.0;
						r_factor_face[j*3+k] = (u[j] - phi_u_star)/(u[neighbor_ids[j*3+k]]-u[j]);
					}
					else
					{
						//cout << "else is here :/\n";
						//r_factor_face[j*3+k] = (2.0 * cell_gradients[j] * barycenter_dist_vector[j*3+k])/(1e-40)-1.0;
						if (fabs(u[j] - phi_u_star)<= 1e-6) r_factor_face[j*3+k] = ((u[j] - phi_u_star) >= 0.0) ? 1.0 : -1.0;
						else r_factor_face[j*3+k] = (u[j] - phi_u_star)/(1e-11);
						//cout << u[j] - phi_u_star << endl;
					}
					//face_values[j*3+k] = u[j] + 0.5 * limiter(r_factor_face[j*3+k], L[j*3+k]) * (u[neighbor_ids[j*3+k]]-u[j]);
					//face_values[j*3+k] = u[j] + limiter(r_factor_face[j*3+k], L[j*3+k])/L[j*3+k] * (u[neighbor_ids[j*3+k]] - u[j]);
					face_values[j*3+k] = u[j] + limiter(r_factor_face[j*3+k], L[j*3+k], u[neighbor_ids[j*3+k]], u[j], (u[j]*(1.0-1.0/L[j*3+k]) + u[neighbor_ids[j*3+k]]*1.0/L[j*3+k])*(face_normals[j*3+k]), skewness_corrector[j*3+k])/L[j*3+k] * (u[neighbor_ids[j*3+k]] - u[j]);
				}
				//else
				//{
				//	face_values[j*3+k] = 0.0;
				//	if (t == 0) cout << "here" << endl;
				//}
			}
			//if ((j == 1570 || j == 1615) && t == 0)
			//{
			//	cout << j << " " << r_factor_face[j*3] << " " << r_factor_face[j*3+1] << " " << r_factor_face[j*3+2] << " " << L[j*3] << " " << L[j*3+1] << " " << L[j*3+2] << endl;
			//}
		}
	}
	for (int j = 0; j < Xrange; j++)
	{
		double temp = 0.0;
		if (!is_tr_boundary[j])
		{
			for (int k = 0; k < 3; k++)
			{
				//if (t==0 && (j==439 || j==568 || j==806 || j==3543))
				//{
				//	cout << j << " " << face_values[k+j*3] << " " << face_values[neighbor_ids[j*3+k]*3] <<  " " << face_values[neighbor_ids[j*3+k]*3+1] <<  " " << face_values[neighbor_ids[j*3+k]*3+2] << endl;
				//}
				if ((face_normals[j*3+k]*velocities[j])*(barycenter_dist_vector[j*3+k]*velocities[j])>=0.0) temp += face_values[k+j*3] * velocities[j] * face_normals[j*3+k]; 
			}
			u[j] = u[j] - dt / tr_areas[j] * temp;
			if (u[j] < 0.0) u[j] = 0.0;
			if (u[j] > 1.0) u[j] = 1.0;
		}
	}
	delete [] face_values;
	delete [] cell_gradients;
	delete [] r_factor_face;
}

void advection_include_2d::calc_is_tr_boundary(vector<Triangle2 *> &triangles)
{
	ofstream data;
	data.open("is_boundary.txt");
	int k = 0;
	for (auto it(triangles.begin()); it != triangles.end(); it++)
	{		
		int how_many = 0;
		Point2 * point;
		for (int i = 0; i < 3; i++)
		{
			point = (*it)->getCorner(i);
			//if (point->x() == 0.0 || point->y() == 0.0 || point->x() == 1.0 || point->y() == 1.0) how_many++;
			if (point->x() == 0.0 || point->y() == 0.0 || point->x() == 1.0 || point->y() == 1.0 || ((point->x()-0.5)*(point->x()-0.5) + (point->y()-0.2)*(point->y()-0.2) < 0.08*0.08) || ((point->x()-0.5)*(point->x()-0.5) + (point->y()-0.5)*(point->y()-0.5) < 0.08*0.08) || ((point->x()-0.5)*(point->x()-0.5) + (point->y()-0.8)*(point->y()-0.8) < 0.08*0.08)) how_many++;
		}
		if (how_many >= 1) is_tr_boundary[k] = true;
		else is_tr_boundary[k] = false;
		data << k << " " << is_tr_boundary[k] << endl;
		k++;
	}
	data.close();
}

void advection_include_2d::calc_vel(vector<Triangle2 *> &triangles)
{
	int k = 0;
	for (auto it(triangles.begin()); it != triangles.end(); it++)
	{
		Point2 point;
		point = (*it)->getBarycenter();
		velocities[k] = Vector2(1.0, 0.0);
		k++;
	}
}

void advection_include_2d::calc_face_normals_and_more(vector<Triangle2 *> &triangles)
{
	ofstream data;
	data.open("face_normals.txt");
	int k = 0;
	double length;
	double x1;
	double y1;
	for (auto it(triangles.begin()); it != triangles.end(); it++)
	{
		if (!is_tr_boundary[k])
		{
			Point2 barry = (*it)->getBarycenter();
			Point2 * corner1;
			Point2 * corner2;
			for (int i = 0; i < 3; i++)
			{
				corner1 = (*it)->getCorner((i+1)%3);
				corner2 = (*it)->getCorner((i+2)%3);
				Point2 center = Point2((corner1->x() + corner2->x())/2, (corner1->y() + corner2->y())/2);
				length = pow(pow(corner1->x()-corner2->x(), 2) + pow(corner1->y() - corner2->y(), 2), 0.5);
				x1 = center.x() + (abs(corner1->y()-corner2->y())*length/pow(pow(corner1->y()-corner2->y(), 2)+pow(corner1->x()-corner2->x(), 2), 0.5)) * ((barry.x() < center.x()) ? 1.0 : -1.0);
				y1 = center.y() + (abs(corner1->x()-corner2->x())*length/pow(pow(corner1->x()-corner2->x(), 2)+pow(corner1->y()-corner2->y(), 2), 0.5)) * ((barry.y() < center.y()) ? 1.0 : -1.0);
				Point2 dest = Point2(x1, y1);

				assert(fabs(length - pow(pow(dest.x()-center.x(), 2) + pow(dest.y()-center.y(),2),0.5))<=1e-10);

				face_normals[k*3+i] = Vector2(dest.x()-center.x(), dest.y()-center.y());
				Point2 neighborbary = (*it)->getOppositeTriangle(i)->getBarycenter();
				L[k*3+i] = (pow(pow(barry.x()-center.x(),2)+pow(barry.y()-center.y(),2),0.5) + pow(pow(neighborbary.x()-center.x(),2)+pow(neighborbary.y()-center.y(),2),0.5))/pow(pow(barry.x()-center.x(),2)+pow(barry.y()-center.y(),2),0.5);
				double y, x;
				if (barry.x() == neighborbary.x())
				{
					x = barry.x();
					y = (x - corner1->x())/(corner2->x()-corner1->x())*(corner2->y()-corner1->y())+corner1->y();
				}
				else if (barry.y() == neighborbary.y())
				{
					y = barry.y();
					x = (y - corner1->y())/(corner2->y()-corner1->y())*(corner2->x()-corner1->x())+corner1->x();
				}
				else if (corner1->x() == corner2->x())
				{
					x = corner1->x();
					y = (x - barry.x())/(neighborbary.x()-barry.x())*(neighborbary.y()-barry.y())+barry.y();
				}
				else if (corner1->y() == corner2->y())
				{
					y = corner1->y();
					x = (y - barry.y())/(neighborbary.y()-barry.y())*(neighborbary.x()-barry.x())+barry.x();
				}
				else
				{ 
					y = (-barry.y()*(corner2->y()-corner1->y())*(neighborbary.x()-barry.x())/((neighborbary.y()-barry.y())*(corner2->x()-corner1->x()))+corner1->y()+(corner2->y()-corner1->y())*(barry.x()-corner1->x())/(corner2->x()-corner1->x()))/(1.0-(corner2->y()-corner1->y())*(neighborbary.x()-barry.x())/((neighborbary.y()-barry.y())*(corner2->x()-corner1->x())));
					x = (neighborbary.x()-barry.x())*(y-barry.y())+barry.x()*(neighborbary.y()-barry.y())/(neighborbary.y()-barry.y());
				}
				skewness_corrector[k*3+i] = Vector2(center.x()-x, center.y()-y);
				//cout << barry.x() << " " << barry.y() << " " << neighborbary.x() << " " << neighborbary.y() << " " << center.x() << " " << center.y() << " " << x << " " << y << endl;
				//cout << barry.x() << " " << barry.y() << " " << center.x() << " " << center.y() << " " << dest.x() << " " << dest.y() << " "<< length << endl;
			}
			data << k << " " << face_normals[k*3].x() << "/" << face_normals[k*3].y() << " " << face_normals[k*3+1].x() << "/" << face_normals[k*3+1].y() << " " << face_normals[k*3+2].x() << "/" << face_normals[k*3+2].y() << endl;
		}
		k++;
	}
	data.close();
}

void advection_include_2d::calc_bary_dist_vec(vector<Triangle2 *> &triangles)
{
	ofstream data;
	data.open("bary_dist.txt");
	int k = 0;
	for (auto it(triangles.begin()); it != triangles.end(); it++)
	{
		if (!is_tr_boundary[k])
		{
			Point2 origin = (*it)->getBarycenter();
			Point2 destination;
			for (int i = 0; i < 3; i++)
			{
				destination = (*it)->getOppositeTriangle(i)->getBarycenter();
				barycenter_dist_vector[k*3+i] = Vector2(destination.x()-origin.x(), destination.y()-origin.y());;
			}
			data << k << " " << barycenter_dist_vector[k*3].x() << "/" << barycenter_dist_vector[k*3].y() << " " << barycenter_dist_vector[k*3+1].x() << "/" << barycenter_dist_vector[k*3+1].y() << " " << barycenter_dist_vector[k*3+2].x() << "/" << barycenter_dist_vector[k*3+2].y() << endl;
		}
		//else for (int i = 0; i < 3; i++) barycenter_dist_vector[k*3+i] = Vector2(0.0, 0.0);
		k++;
	}
	data.close();
}

void advection_include_2d::calc_neighbor_ids(vector<Triangle2 *> &triangles)
{
	ofstream data;
	data.open("neighbor_ids.txt");
	int k = 0;
	for (auto it(triangles.begin()); it != triangles.end(); it++)
	{
		if (!is_tr_boundary[k])
		{
			for (int i = 0; i < 3; i++)
			{
				auto neighbor = find(triangles.begin(), triangles.end(),(*it)->getOppositeTriangle(i));
				if (neighbor != triangles.end()) neighbor_ids[k*3+i] = distance(triangles.begin(), neighbor);
				else cout << "something went terribly wrong" << endl;
			}
		}
		else for (int i = 0; i < 3; i++) neighbor_ids[k*3+i] = -1;
		data << k << " " << neighbor_ids[k*3] << " " << neighbor_ids[k*3+1] << " " << neighbor_ids[k*3+2] << endl;
		k++;
	}
	//for (int i = 0; i < Xrange; i++)
	//{
		//cout << i << " " << neighbor_ids[i*3] << " " << neighbor_ids[i*3+1] << " " << neighbor_ids[i*3+2] << " " << neighbor_ids[neighbor_ids[i*3]*3] << " " << neighbor_ids[neighbor_ids[i*3]*3+1]<< " " << neighbor_ids[neighbor_ids[i*3]*3+2] << endl;
	//}
	data.close();
}

void advection_include_2d::calc_tr_areas(vector<Triangle2 *> &triangles)
{
	ofstream data;
	data.open("tr_areas.txt");
	double area_min = 100.0;
	double denom_max = 0.0;
	int k = 0;
	for (auto it(triangles.begin()); it != triangles.end(); it++)
	{
		tr_areas[k] = (*it)->getArea2D();
		if (area_min > tr_areas[k]) area_min = tr_areas[k];
		for (int i = 0; i < 3; i++) if (velocities[k]*face_normals[k*3+i] > denom_max) denom_max = velocities[k]*face_normals[k*3+i];
		data << k << " " << tr_areas[k] << endl;
		k++;
	}
	cout << "recommended timestep " << area_min/denom_max << " or less" << endl;
	data.close();
	
}

double advection_include_2d::limiter(double r_factor, double L, double value_downwind, double value_center, Vector2 face_gradient, Vector2 distance)
{
	//return 1.0;
	//double result;
	double temp_result;
	temp_result = (0.5*L*r_factor+0.5*L*fabs(r_factor))/(L-1+fabs(r_factor));//van leer
	//temp_result = max(0.0, min(1.0, r_factor));//minmod
	//return temp_result;//minmod
	//temp_result = max(min(1.0, L*r_factor), min(r_factor, L));//superbee non custom
	//result = max(temp_result, 0.0) + L*face_gradient*distance/(value_downwind - value_center);
	//result = temp_result + L*face_gradient*distance/(value_downwind - value_center);//van leer	
	//temp_result = max(min(1.0, L*result), min(result, L));
	//if (result < 0.0) result = 0.0;
	//else if (result > max(0.0, min(L*r_factor, L))) result = max(0.0, min(L*r_factor, L));
	//else result = max(0.0, min(L*r_factor, L));	
	//return max(temp_result, 0.0);//superbee
	//return result;
	return temp_result;//van leer
}
class CustomParameters:public MeshGenParams
{
public:
        CustomParameters(Zone2* pZone):MeshGenParams(pZone)
        {
        }
        double getMaxEdgeLength(Triangle2* pT)
        {
                Point2 barycenter(pT->getBarycenter());
                if((pow(barycenter.x()-0.5, 2) + pow(barycenter.y()-0.5, 2) <= pow(0.2, 2)) || (barycenter.y()>0.95 || barycenter.x()>0.95 || barycenter.y()<0.05 || barycenter.x()<0.05))
                {
                        return 0.05;
                }
                else
                {
                        return 0.15;
                }
        }
};

int main()
{
	Fade_2D dt;
	vector<Point2> square{Point2(0.0, 0.0), Point2(0.0, 1.0), Point2(1.0, 1.0), Point2(1.0, 0.0)};
	dt.insert(square);
	//vector<Point2> circle;
	vector<Point2> circle1;
	vector<Point2> circle2;
	vector<Point2> circle3;
	generateCircle(10, 0.5, 0.2, 0.075, 0.075, circle1);
	generateCircle(10, 0.5, 0.5, 0.075, 0.075, circle2);
	generateCircle(10, 0.5, 0.8, 0.075, 0.075, circle3);
	//generateCircle(3, 0.5, 0.5, 0.01, 0.01, circle);
	//generateCircle(25, 0.5, 0.5, 0.1, 0.1, circle);
	//vector<Segment2> circle_segment;
	//vector<Segment2> square_segment;
	vector<Segment2> circle_segment1;
	vector<Segment2> circle_segment2;
	vector<Segment2> circle_segment3;
	//for (int i = 0; i < circle.size() - 1; i++)
	//{
	//	circle_segment.push_back(Segment2(circle[i], circle[i+1]));
	//}
	for (int i = 0; i < circle1.size() - 1; i++)
	{
		circle_segment1.push_back(Segment2(circle1[i], circle1[i+1]));
		circle_segment2.push_back(Segment2(circle2[i], circle2[i+1]));
		circle_segment3.push_back(Segment2(circle3[i], circle3[i+1]));
	}
	//for (int i = 0; i < square.size() - 1; i++)
	//{
	//	square_segment.push_back(Segment2(square[i], square[i+1]));
	//}
	//circle_segment.push_back(Segment2(circle.back(), circle.front()));
	circle_segment1.push_back(Segment2(circle1.back(), circle1.front()));
	circle_segment2.push_back(Segment2(circle2.back(), circle2.front()));
	circle_segment3.push_back(Segment2(circle3.back(), circle3.front()));
	//square_segment.push_back(Segment2(square.back(), square.front()));
	//ConstraintGraph2 * pCG = dt.createConstraint(circle_segment, CIS_CONSTRAINED_DELAUNAY);
	ConstraintGraph2 * pCG1 = dt.createConstraint(circle_segment1, CIS_CONSTRAINED_DELAUNAY);
	ConstraintGraph2 * pCG2 = dt.createConstraint(circle_segment2, CIS_CONSTRAINED_DELAUNAY);
	ConstraintGraph2 * pCG3 = dt.createConstraint(circle_segment3, CIS_CONSTRAINED_DELAUNAY);
	//ConstraintGraph2 * pCG = dt.createConstraint(square_segment, CIS_CONSTRAINED_DELAUNAY);	
	dt.applyConstraintsAndZones();
	//Zone2 * pZoneOutside(dt.createZone(pCG, ZL_OUTSIDE));
	//Zone2 * pZoneOutside(dt.createZone(pCG, ZL_INSIDE));
	//Zone2 * pZoneOutside(dt.createZone(NULL,ZL_GLOBAL));
	Zone2 * pZone1(dt.createZone(pCG1, ZL_OUTSIDE));
	Zone2 * pZone2(dt.createZone(pCG2, ZL_OUTSIDE));
	Zone2 * pZone3(dt.createZone(pCG3, ZL_OUTSIDE));
	Zone2 * pZoneTemp(zoneIntersection(pZone1,pZone2));
	Zone2 * pZoneOutside(zoneIntersection(pZone3,pZoneTemp));
	dt.refine(pZoneOutside->convertToBoundedZone(), 30, 0.01, 0.03, true);
	//CustomParameters params(pZoneOutside->convertToBoundedZone());
	//MeshGenParams params(pZoneOutside->convertToBoundedZone());
	//params.minAngleDegree = 30;
	//params.minEdgeLength = 0.1;
	//params.maxEdgeLength = 0.15;
	//params.growFactor = 3;
	//dt.refineAdvanced(&params);

	Visualizer2 visualizer("zone_visualizer.ps");
	vector<Triangle2*> visual;
	pZoneOutside->getTriangles(visual);
	for (auto it(visual.begin()); it != visual.end(); it++)
	{
		visualizer.addObject(**it, Color(1,0,0,0.02, true));
	}
	//pZoneOutside->show("zone.ps", true, true);
	//dt.show("square.ps");
	visualizer.writeFile();
	//vector<Point2*> vertices;
	//pZoneOutside->convertToBoundedZone()->getVertices(vertices);
	cout << "THE NUMBER OF triangles IS " << visual.size() << endl;
	unique_ptr<advection_include_2d> equation = make_unique<advection_include_2d>(visual, 0.0001);
	ofstream data;
	data.open("data.txt");
	data << fixed << setprecision(6);
	double error = 0.0;
	for (int t = 0; t < 10000; t++)
	{
		int k = 0;
		
		for (auto it(visual.begin()); it != visual.end(); it++)
		{
			double diff = 0.0;
			if ((t<10000) &&((*it)->getBarycenter().x() < 0.9) && ((*it)->getBarycenter().y() < 0.9) && ((*it)->getBarycenter().y() > 0.1) && ((*it)->getBarycenter().x() > 0.1)) diff = fabs(equation->u[k] - exp(-pow((*it)->getBarycenter().y()-0.5,2)*200.0-pow((*it)->getBarycenter().x()-0.2 - t*0.0001,2)*200.0));
			if (diff>error) error = diff;
			if (t%100 == 0) data <<(*it)->getBarycenter().x() << " " << (*it)->getBarycenter().y() << " " << equation->u[k] << " " << k << endl;
			//if (t%20 == 0) data <<(*it)->getBarycenter().x() << " " << (*it)->getBarycenter().y() << " " << exp(-pow((*it)->getBarycenter().y()-0.5,2)*200.0-pow((*it)->getBarycenter().x()-0.2 - t*0.0005,2)*200.0) <<endl;			
			k++;
		}
		if (t%100 == 0) data << endl << endl;
		equation->solver(t);
	}
	cout << error << " is the error\n";
	data.close();
	return 0;
}

