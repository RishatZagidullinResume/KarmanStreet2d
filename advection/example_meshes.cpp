#include "example_meshes.h"

namespace example_meshes
{
	void mesh_three_circles(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length)
	{
		vector<Point2> square{Point2(0.0, 0.0), Point2(0.0, 1.0), Point2(1.0, 1.0), Point2(1.0, 0.0)};
		dt.insert(square);
		vector<Point2> circle1;
		vector<Point2> circle2;
		vector<Point2> circle3;
		generateCircle(10, 0.5, 0.2, 0.075, 0.075, circle1);
		generateCircle(10, 0.5, 0.5, 0.075, 0.075, circle2);
		generateCircle(10, 0.5, 0.8, 0.075, 0.075, circle3);
		vector<Segment2> square_segment;
		vector<Segment2> circle_segment1;
		vector<Segment2> circle_segment2;
		vector<Segment2> circle_segment3;
		for (int i = 0; i < circle1.size() - 1; i++)
		{
			circle_segment1.push_back(Segment2(circle1[i], circle1[i+1]));
			circle_segment2.push_back(Segment2(circle2[i], circle2[i+1]));
			circle_segment3.push_back(Segment2(circle3[i], circle3[i+1]));
		}
		for (int i = 0; i < square.size() - 1; i++)
		{
			square_segment.push_back(Segment2(square[i], square[i+1]));
		}
		circle_segment1.push_back(Segment2(circle1.back(), circle1.front()));
		circle_segment2.push_back(Segment2(circle2.back(), circle2.front()));
		circle_segment3.push_back(Segment2(circle3.back(), circle3.front()));
		square_segment.push_back(Segment2(square.back(), square.front()));
		ConstraintGraph2 * pCG1 = dt.createConstraint(circle_segment1, CIS_CONSTRAINED_DELAUNAY);
		ConstraintGraph2 * pCG2 = dt.createConstraint(circle_segment2, CIS_CONSTRAINED_DELAUNAY);
		ConstraintGraph2 * pCG3 = dt.createConstraint(circle_segment3, CIS_CONSTRAINED_DELAUNAY);
		ConstraintGraph2 * pCG = dt.createConstraint(square_segment, CIS_CONSTRAINED_DELAUNAY);	
		dt.applyConstraintsAndZones();
		Zone2 * pZone1(dt.createZone(pCG1, ZL_OUTSIDE));
		Zone2 * pZone2(dt.createZone(pCG2, ZL_OUTSIDE));
		Zone2 * pZone3(dt.createZone(pCG3, ZL_OUTSIDE));
		Zone2 * pZoneTemp(zoneIntersection(pZone1,pZone2));
		Zone2 * pZoneOutside(zoneIntersection(pZone3,pZoneTemp));
		dt.refine(pZoneOutside->convertToBoundedZone(), 30, min_length, max_length, true);
		pZoneOutside->getTriangles(triangles);
		pZoneOutside->convertToBoundedZone()->getVertices(nodes);	
	}

	void mesh_one_circle(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length)
	{
		vector<Point2> square{Point2(0.0, 0.0), Point2(0.0, 0.5), Point2(1.5, 0.5), Point2(1.5, 0.0)};
		dt.insert(square);
		vector<Point2> circle;
		generateCircle(20, 0.3, 0.25, 0.04, 0.06, circle);
		vector<Segment2> circle_segment;
		vector<Segment2> square_segment;
		for (int i = 0; i < circle.size() - 1; i++)
		{
			circle_segment.push_back(Segment2(circle[i], circle[i+1]));
		}
		for (int i = 0; i < square.size() - 1; i++)
		{
			square_segment.push_back(Segment2(square[i], square[i+1]));
		}
		circle_segment.push_back(Segment2(circle.back(), circle.front()));
		square_segment.push_back(Segment2(square.back(), square.front()));
		ConstraintGraph2 * pCG = dt.createConstraint(circle_segment, CIS_CONSTRAINED_DELAUNAY);
		dt.applyConstraintsAndZones();
		Zone2 * pZoneOutside(dt.createZone(pCG, ZL_OUTSIDE));
		dt.refine(pZoneOutside->convertToBoundedZone(), 30, min_length, max_length, true);
		pZoneOutside->getTriangles(triangles);
		pZoneOutside->convertToBoundedZone()->getVertices(nodes);
	}

	void mesh_unstructured_square(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length)
	{
		vector<Point2> square{Point2(0.0, 0.0), Point2(0.0, 1.0), Point2(1.0, 1.0), Point2(1.0, 0.0)};
		dt.insert(square);
		vector<Point2> circle;
		generateCircle(3, 0.5, 0.5, 0.01, 0.01, circle);
		vector<Segment2> circle_segment;
		vector<Segment2> square_segment;
		for (int i = 0; i < circle.size() - 1; i++)
		{
			circle_segment.push_back(Segment2(circle[i], circle[i+1]));
		}
		for (int i = 0; i < square.size() - 1; i++)
		{
			square_segment.push_back(Segment2(square[i], square[i+1]));
		}
		circle_segment.push_back(Segment2(circle.back(), circle.front()));
		square_segment.push_back(Segment2(square.back(), square.front()));
		ConstraintGraph2 * pCG = dt.createConstraint(circle_segment, CIS_CONSTRAINED_DELAUNAY);
		dt.applyConstraintsAndZones();
		Zone2 * pZoneOutside(dt.createZone(NULL,ZL_GLOBAL));
		dt.refine(pZoneOutside->convertToBoundedZone(), 30, min_length, max_length, true);
		pZoneOutside->getTriangles(triangles);
		pZoneOutside->convertToBoundedZone()->getVertices(nodes);
	
	}

	void mesh_unstructured_rect(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length)
	{
		vector<Point2> square{Point2(0.0, 0.0), Point2(0.0, 0.5), Point2(2.0, 0.5), Point2(2.0, 0.0)};
		dt.insert(square);
		vector<Point2> circle;
		generateCircle(6, 1.0, 0.25, 0.025, 0.025, circle);
		vector<Segment2> circle_segment;
		vector<Segment2> square_segment;
		for (int i = 0; i < circle.size() - 1; i++)
		{
			circle_segment.push_back(Segment2(circle[i], circle[i+1]));
		}
		for (int i = 0; i < square.size() - 1; i++)
		{
			square_segment.push_back(Segment2(square[i], square[i+1]));
		}
		circle_segment.push_back(Segment2(circle.back(), circle.front()));
		square_segment.push_back(Segment2(square.back(), square.front()));
		ConstraintGraph2 * pCG = dt.createConstraint(circle_segment, CIS_CONSTRAINED_DELAUNAY);
		dt.applyConstraintsAndZones();
		Zone2 * pZoneOutside(dt.createZone(NULL,ZL_GLOBAL));
		dt.refine(pZoneOutside->convertToBoundedZone(), 30, min_length, max_length, true);
		pZoneOutside->getTriangles(triangles);
		pZoneOutside->convertToBoundedZone()->getVertices(nodes);
	
	}

	void mesh_structured_rect_var(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length, double x_length, double y_length)
	{
		vector<Point2> square{Point2(0.0, 0.0), Point2(0.0, y_length), Point2(x_length, y_length), Point2(x_length, 0.0)};
		dt.insert(square);
		vector<Segment2> square_segment;
		for (int i = 0; i < square.size() - 1; i++)
		{
			square_segment.push_back(Segment2(square[i], square[i+1]));
		}
		square_segment.push_back(Segment2(square.back(), square.front()));
		ConstraintGraph2 * pCG = dt.createConstraint(square_segment, CIS_CONSTRAINED_DELAUNAY);	
		dt.applyConstraintsAndZones();
		Zone2 * pZoneOutside(dt.createZone(pCG, ZL_INSIDE));
		dt.refine(pZoneOutside->convertToBoundedZone(), 30, min_length, max_length, true);
		pZoneOutside->getTriangles(triangles);
		pZoneOutside->convertToBoundedZone()->getVertices(nodes);
	}

	void mesh_little_unstructured(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length)
	{
		vector<Point2> square{Point2(0.0, 0.0), Point2(0.0, 0.2), Point2(0.2, 0.2), Point2(0.2, 0.0)};
		dt.insert(square);
		vector<Point2> circle;
		generateCircle(3, 0.05, 0.05, 0.01, 0.01, circle);
		vector<Segment2> circle_segment;
		vector<Segment2> square_segment;
		for (int i = 0; i < circle.size() - 1; i++)
		{
			circle_segment.push_back(Segment2(circle[i], circle[i+1]));
		}
		for (int i = 0; i < square.size() - 1; i++)
		{
			square_segment.push_back(Segment2(square[i], square[i+1]));
		}
		circle_segment.push_back(Segment2(circle.back(), circle.front()));
		square_segment.push_back(Segment2(square.back(), square.front()));
		ConstraintGraph2 * pCG = dt.createConstraint(circle_segment, CIS_CONSTRAINED_DELAUNAY);
		dt.applyConstraintsAndZones();
		Zone2 * pZoneOutside(dt.createZone(NULL,ZL_GLOBAL));
		dt.refine(pZoneOutside->convertToBoundedZone(), 30, min_length, max_length, true);
		pZoneOutside->getTriangles(triangles);
		pZoneOutside->convertToBoundedZone()->getVertices(nodes);
	}

	void mesh_tube(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length)
	{
		vector<Point2> square{Point2(0.0, 0.0), Point2(0.2, 0.0), Point2(0.4, 0.1), Point2(1.0, 0.1), Point2(1.0, 0.3), Point2(0.4, 0.3), Point2(0.2, 0.2), Point2(0.0, 0.2)};
		dt.insert(square);
		vector<Segment2> square_segment;
		for (int i = 0; i < square.size() - 1; i++)
		{
			square_segment.push_back(Segment2(square[i], square[i+1]));
		}
		square_segment.push_back(Segment2(square.back(), square.front()));
		ConstraintGraph2 * pCG = dt.createConstraint(square_segment, CIS_CONSTRAINED_DELAUNAY);	
		dt.applyConstraintsAndZones();
		Zone2 * pZoneOutside(dt.createZone(pCG, ZL_INSIDE));
		dt.refine(pZoneOutside->convertToBoundedZone(), 30, min_length, max_length, true);
		pZoneOutside->getTriangles(triangles);
		pZoneOutside->convertToBoundedZone()->getVertices(nodes);
	}

	void setuha(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length, int parabola_points_num)
	{
		vector<Point2> square{Point2(0.0, 0.0), Point2(0.0, 0.5)};
		double dx = 1.0/parabola_points_num;
		for (int i = 0; i <= parabola_points_num; i++)
		{
			double x = 2.0+dx*i;
			if (i == parabola_points_num && x != 3.0)
			{
				cout << "x forced 3.0" << endl;
				x = 3.0;
			}
			square.push_back(Point2(x, -0.125*x*x+0.5*x));
		}
		for (int i = parabola_points_num; i >= 0; i--)
		{
			double x = 2.0+dx*i;
			if (i == parabola_points_num && x != 3.0)
			{
				cout << "x forced 3.0" << endl;
				x = 3.0;
			}
			square.push_back(Point2(x, 0.125*x*x-0.5*x+0.5));
		}
		dt.insert(square);
		vector<Segment2> square_segment;
		for (int i = 0; i < square.size() - 1; i++)
		{
			square_segment.push_back(Segment2(square[i], square[i+1]));
		}
		square_segment.push_back(Segment2(square.back(), square.front()));
		ConstraintGraph2 * pCG = dt.createConstraint(square_segment, CIS_CONSTRAINED_DELAUNAY);	
		dt.applyConstraintsAndZones();
		Zone2 * pZoneOutside(dt.createZone(pCG, ZL_INSIDE));
		dt.refine(pZoneOutside->convertToBoundedZone(), 30, min_length, max_length, true);
		pZoneOutside->getTriangles(triangles);
		pZoneOutside->convertToBoundedZone()->getVertices(nodes);
	}

	vector<double> mesh_structured_rect_var_mpi(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length, double x_length, double y_length, int n_processes)
	{
		vector<Point2> square{Point2(0.0, 0.0), Point2(0.0, y_length), Point2(x_length, y_length), Point2(x_length, 0.0)};
		dt.insert(square);
		vector<Segment2> square_segment;
		for (int i = 0; i < square.size() - 1; i++)
		{
			square_segment.push_back(Segment2(square[i], square[i+1]));
		}
		square_segment.push_back(Segment2(square.back(), square.front()));
		ConstraintGraph2 * pCG = dt.createConstraint(square_segment, CIS_CONSTRAINED_DELAUNAY);	

		double x_len = x_length;
		double y_len = y_length;
		while ( (y_length/y_len)*(x_length/x_len) != n_processes)
		{
			
			x_len > y_len ? x_len/=2.0 : y_len/=2.0;
		}
		vector<ConstraintGraph2 *> procCG;
		for (int i = 0; i < x_length/x_len - 1; i++)
		{
			Point2 point1 = Point2((i+1)*x_len, 0);
			Point2 point2 = Point2((i+1)*x_len, y_length);
			vector<Segment2> segment{Segment2(point1, point2)};
			procCG.push_back(dt.createConstraint(segment, CIS_CONFORMING_DELAUNAY));	
		}
		for (int i = 0; i < y_length/y_len - 1; i++)
		{
			Point2 point1 = Point2(0, (i+1)*y_len);
			Point2 point2 = Point2(x_length, (i+1)*y_len);
			vector<Segment2> segment{Segment2(point1, point2)};
			procCG.push_back(dt.createConstraint(segment, CIS_CONFORMING_DELAUNAY));
		}
		dt.applyConstraintsAndZones();
		Zone2 * pZoneOutside(dt.createZone(pCG, ZL_INSIDE));
		dt.refine(pZoneOutside->convertToBoundedZone(), 30, min_length, max_length, true);
		pZoneOutside->getTriangles(triangles);
		pZoneOutside->convertToBoundedZone()->getVertices(nodes);
		return {x_len, y_len};
	}

	void setuha_mpi(Fade_2D& dt, vector<Triangle2*>& triangles, vector<Point2*>& nodes, double min_length, double max_length, int parabola_points_num, int n_processes)
	{
		vector<Point2> square{Point2(0.0, 0.0), Point2(0.0, 0.5)};
		double dx = 1.0/parabola_points_num;
		for (int i = 0; i <= parabola_points_num; i++)
		{
			double x = 2.0+dx*i;
			if (i == parabola_points_num && x != 3.0)
			{
				cout << "x forced 3.0" << endl;
				x = 3.0;
			}
			square.push_back(Point2(x, -0.125*x*x+0.5*x));
		}
		for (int i = parabola_points_num; i >= 0; i--)
		{
			double x = 2.0+dx*i;
			if (i == parabola_points_num && x != 3.0)
			{
				cout << "x forced 3.0" << endl;
				x = 3.0;
			}
			square.push_back(Point2(x, 0.125*x*x-0.5*x+0.5));
		}
		dt.insert(square);
		vector<Segment2> square_segment;
		for (int i = 0; i < square.size() - 1; i++)
		{
			square_segment.push_back(Segment2(square[i], square[i+1]));
		}
		square_segment.push_back(Segment2(square.back(), square.front()));
		ConstraintGraph2 * pCG = dt.createConstraint(square_segment, CIS_CONSTRAINED_DELAUNAY);	
		
		double x_length = 3.0;
		double y_length = 0.5;
		double x_len = 3.0;
		double y_len = 0.5;
		while ( (y_length/y_len)*(x_length/x_len) != n_processes)
		{
			
			x_len > y_len ? x_len/=2.0 : y_len/=2.0;
		}
		vector<ConstraintGraph2 *> procCG;
		for (int i = 0; i < x_length/x_len - 1; i++)
		{
			Point2 point1 = Point2((i+1)*x_len, 0);
			Point2 point2 = Point2((i+1)*x_len, y_length);
			vector<Segment2> segment{Segment2(point1, point2)};
			procCG.push_back(dt.createConstraint(segment, CIS_CONFORMING_DELAUNAY));	
		}
		for (int i = 0; i < y_length/y_len - 1; i++)
		{
			Point2 point1 = Point2(0, (i+1)*y_len);
			Point2 point2 = Point2(x_length, (i+1)*y_len);
			vector<Segment2> segment{Segment2(point1, point2)};
			procCG.push_back(dt.createConstraint(segment, CIS_CONFORMING_DELAUNAY));
		}

		dt.applyConstraintsAndZones();
		Zone2 * pZoneOutside(dt.createZone(pCG, ZL_INSIDE));
		dt.refine(pZoneOutside->convertToBoundedZone(), 30, min_length, max_length, true);
		pZoneOutside->getTriangles(triangles);
		pZoneOutside->convertToBoundedZone()->getVertices(nodes);
	}
}


