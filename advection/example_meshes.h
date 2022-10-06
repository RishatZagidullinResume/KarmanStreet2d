#include "../fade2d/Fade_2D.h"

using namespace GEOM_FADE2D;
using namespace std;

namespace example_meshes
{
	void mesh_three_circles(Fade_2D&, vector<Triangle2*>&, vector<Point2*>&, double, double);
	void mesh_one_circle(Fade_2D&, vector<Triangle2*>&, vector<Point2*>&, double, double);
	void mesh_unstructured_square(Fade_2D&, vector<Triangle2*>&, vector<Point2*>&, double, double);
	void mesh_unstructured_rect(Fade_2D&, vector<Triangle2*>&, vector<Point2*>&, double, double);
	void mesh_tube(Fade_2D&, vector<Triangle2*>&, vector<Point2*>&, double, double);
	void mesh_little_unstructured(Fade_2D&, vector<Triangle2*>&, vector<Point2*>&, double, double);
	void setuha(Fade_2D&, vector<Triangle2*>&, vector<Point2*>&, double, double, int parabola_points_num);
	void mesh_structured_rect_var(Fade_2D&, vector<Triangle2*>&, vector<Point2*>&, double, double, double, double);
	vector<double> mesh_structured_rect_var_mpi(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length, double x_length, double y_length, int n_processes);
	void setuha_mpi(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length, int parabola_points_num, int n_processes);

}
