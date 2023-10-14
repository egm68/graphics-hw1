////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;

double inline det(const Point &u, const Point &v, const Point &a, const Point &b) {
	// TODO
    double determinant = (imag(v) - imag(u)) * (real(a) - real(b)) - (imag(b) - imag(a)) * (real(u) - real(v));

	return determinant;
}

// Return true iff [a,b] intersects [c,d], and store the intersection in ans
bool intersect_segment(const Point &a, const Point &b, const Point &c, const Point &d, Point &ans) {
	std::vector<Point> seg_points;
	seg_points.push_back(a);
	seg_points.push_back(b);
	seg_points.push_back(c);
	seg_points.push_back(d);

	// TODO
	double determinant = det(a, b, c, d);
	if (determinant == 0)
    {
	   std::cout << "DET 0" << "\n";
       return false;
    }
    else
    {
		double intersection_x = ((real(c) - real(d)) * ((imag(b) - imag(a)) * (real(a)) + (real(a) - real(b)) * (imag(a))) - (real(a) - real(b)) * ((imag(d) - imag(c)) * (real(c)) + (real(c) - real(d)) * (imag(c)))) / determinant;
        double intersection_y = ((imag(b) - imag(a)) * ((imag(d) - imag(c)) * (real(c)) + (real(c) - real(d)) * (imag(c))) - (imag(d) - imag(c)) * ((imag(b) - imag(a)) * (real(a)) + (real(a) - real(b)) * (imag(a)))) / determinant;
		
		double maxXsegments = -10000000000000000.0;
		double minXsegments = 10000000000000000.0;
		double maxYsegments = -10000000000000000.0;
		double minYsegments = 10000000000000000.0;
		for (size_t i = 0; i < seg_points.size(); ++i) {
			if (real(seg_points[i]) < minXsegments){
				minXsegments = real(seg_points[i]);
			}
			if (real(seg_points[i]) > maxXsegments){
				maxXsegments = real(seg_points[i]);
			}
			if (imag(seg_points[i]) < minYsegments){
				minYsegments = imag(seg_points[i]);
			}
			if (imag(seg_points[i]) > maxYsegments){
				maxYsegments = imag(seg_points[i]);
			}
		}

		if (intersection_x >= minXsegments && intersection_x <= maxXsegments && intersection_y >= minYsegments && intersection_y <= maxYsegments){
			ans = Point(intersection_x, intersection_y);
        	return true;
		}
		else{
			return false;
		}
    }

}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const Polygon &poly, const Point &query) {
	// 1. Compute bounding box 
	// TODO
	Polygon bounding_box;
	double minX = 10000000000000000.0;
	double maxX = -10000000000000000.0;
	double minY = 10000000000000000.0;
	double maxY = -10000000000000000.0;
	for (size_t i = 0; i < poly.size(); ++i) {
		if (real(poly[i]) < minX){
			minX = real(poly[i]);
		}
		if (real(poly[i]) > maxX){
			maxX = real(poly[i]);
		}
		if (imag(poly[i]) < minY){
			minY = imag(poly[i]);
		}
		if (imag(poly[i]) > maxY){
			maxY = imag(poly[i]);
		}
	}
	//clockwise starting from top left
	bounding_box.insert(bounding_box.end(), Point(minX, maxY)); 
	bounding_box.insert(bounding_box.end(), Point(maxX, maxY)); 
	bounding_box.insert(bounding_box.end(), Point(maxX, minY)); 
	bounding_box.insert(bounding_box.end(), Point(minX, minY)); 

	//set coordinate of a point outside the polygon
	Point outside(5.0*maxX, 5.0*maxY);

	// 2. Cast a ray from the query point to the 'outside' point, count number of intersections
	// TODO
	std::vector<Point> intersections;
	for(int i = 0; i < poly.size() - 1; ++i){ //for each edge of polygon
		Point ans(0.0, 0.0);
		if (intersect_segment(query, outside, poly[i], poly[i+1], ans) == true){
			intersections.push_back(ans);
		}
	}

	std::cout << "#INTERSECTIONS: " << "\n";
	std::cout << intersections.size() << "\n";

	if (intersections.size()%2 == 0){
		return false;
	}
	else{
		return true;
	}
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Point> load_xyz(const std::string &filename) {
	std::vector<Point> points;
	std::ifstream in(filename);
	// TODO
	int counter = 0;
	for( std::string line; getline( in, line ); )
	{
		if(counter != 0){
			int delimCounter = 0;
			Point newReadPoint = Point(0.0, 0.0);
			int spaceIndex = 0;
			for (int i = 0; i < line.size(); i++) {
				if (line.substr(i,1) == " " && delimCounter == 0){
					newReadPoint.real(std::stod(line.substr(0,i)));
					delimCounter = 1;
					spaceIndex = i;
				}
				else if (line.substr(i,1) == " " && delimCounter == 1){
					newReadPoint.imag(std::stod(line.substr(spaceIndex+1,i)));
				}
			}
			points.insert(points.end(),newReadPoint);
		}
		counter = counter + 1;
	}

	std::cout << "LENGTH OF POINTS VECTOR: " << points.size() << "\n";

	std::cout << "FIRST 5 POINTS: " << "\n";
	std::cout << points[0] << "\n";
	std::cout << points[1] << "\n";
	std::cout << points[2] << "\n";
	std::cout << points[3] << "\n";
	std::cout << points[4] << "\n";

	return points;
}

Polygon load_obj(const std::string &filename) {
	std::vector<Point> poly;
	std::ifstream in(filename);
	// TODO
	for( std::string line; getline( in, line ); ){
		if (line.substr(0,1) == "v"){
			int delimCounter = 0;
			Point newReadPoint = Point(0.0, 0.0);
			int spaceIndex = 0;
			for (int i = 2; i < line.size(); i++) {
				if (line.substr(i,1) == " " && delimCounter == 0){
					newReadPoint.real(std::stod(line.substr(2,i)));
					delimCounter = 1;
					spaceIndex = i;
				}
				else if (line.substr(i,1) == " " && delimCounter == 1){
					newReadPoint.imag(std::stod(line.substr(spaceIndex+1,i)));
				}
			}
			poly.insert(poly.end(),newReadPoint);
		}

	}

	std::cout << "LENGTH OF POLY VECTOR: " << poly.size() << "\n";

	std::cout << "FIRST 5 poly: " << "\n";
	std::cout << poly[0] << "\n";
	std::cout << poly[1] << "\n";
	std::cout << poly[2] << "\n";
	std::cout << poly[3] << "\n";
	std::cout << poly[4] << "\n";

	return poly;
}

void save_xyz(const std::string &filename, const std::vector<Point> &result) {
	// TODO
	std::ofstream out(filename);
	if (!out.is_open()) {
		throw std::runtime_error("failed to open file " + filename);
	}
	out << std::fixed;
	out << result.size() << "\n";
	for (const auto &v : result) {
		out << v.real() << ' ' << v.imag() << " 0\n";
	}
	out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
	if (argc <= 3) {
		std::cerr << "Usage: " << argv[0] << " points.xyz poly.obj result.xyz" << std::endl;
	}
	std::cout << "HELLO WORLD PIP" << "\n";
	std::vector<Point> points = load_xyz("/Users/erinmcgowan/Archive/points.xyz");
	Polygon poly = load_obj("/Users/erinmcgowan/Archive/polygon.obj");
	std::vector<Point> result;
	for (size_t i = 0; i < points.size(); ++i) {
		if (is_inside(poly, points[i])) {
			result.push_back(points[i]);
		}
	}
	std::cout << "NUMBER OF POINTS INSIDE: " << "\n";
	std::cout << result.size() << "\n";
	save_xyz(argv[3], result);
	return 0;
}
