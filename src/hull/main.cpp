////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>
#include <sstream>
#include <string>
////////////////////////////////////////////////////////////////////////////////

typedef std::complex<double> Point;
typedef std::vector<Point> Polygon;

//task 1 or 3?
//returns distance between two points
//sqrt((P1x - P2x)2 + (P1y - P2y)2)
double inline det(const Point &u, const Point &v) {
	// TODO
	return sqrt(pow(real(u) - real(v), 2) + pow(imag(u) - imag(v), 2));
}

//task 2
struct Compare {
	Point p0; // Leftmost point of the poly

	bool operator ()(const Point &p1, const Point &p2) {
		// TODO
		double rotatedX1 = imag(p1);
		double rotatedY1 = -1 * real(p1);
		double rotatedX2 = imag(p2);
		double rotatedY2 = -1 * real(p2);
		double rotatedX0 = imag(p0);
		double rotatedY0 = -1 * real(p0); 

		double angle1 = -1 * atan2(rotatedX1 - rotatedX0, -1 * (rotatedY1 - rotatedY0));
		double angle2 = -1 * atan2(rotatedX2 - rotatedX0, -1 * (rotatedY2 - rotatedY0));

		if (angle1 < angle2){
			return true;
		}
		else{
			return false;
		}
	}
};

//task 1
//arccos((P122 + P132 - P232) / (2 * P12 * P13))
bool inline salientAngle(Point &a, Point &b, Point &c) {
	// TODO
    int angle = (imag(b) - imag(a)) * (real(c) - real(b)) - (real(b) - real(a)) * (imag(c) - imag(b));
 
	if (angle > 0){
		return false;
	}
	else{
		return true;
	}
}

////////////////////////////////////////////////////////////////////////////////
//task 3
Polygon convex_hull(std::vector<Point> &points) {
	Compare order;
	// TODO
	//find point(s) with the lowest y-coordinate. 
	std::vector<Point> minYpoints;
	double minY = 100000000000000.0;
	double minX = 100000000000000.0;
	for (int i = 0; i < points.size(); i++) {
  		//cout << i << "\n";
		if (imag(points[i]) < minY){ //new minY
			minYpoints.clear();
			minY = imag(points[i]);
			minYpoints.insert(minYpoints.end(),points[i]);
		}
		else if (imag(points[i]) == minY){ //second minY
			minYpoints.insert(minYpoints.end(),points[i]);
		}
	}
	//If the lowest y-coordinate exists in more than one point in the set, 
	//the point with the lowest x-coordinate out of the candidates should be chosen. 
	std::vector<Point> minXYpoints;
	if (minYpoints.size() > 1){
		for (int i = 0; i < minYpoints.size(); i++) {
			if (real(minYpoints[i]) < minX){
				minXYpoints.clear();
				minX = real(minYpoints[i]);
				minXYpoints.insert(minXYpoints.end(),minYpoints[i]);
			}
		}
	}
	else{
		minXYpoints = minYpoints;
	}
	order.p0 = minXYpoints[0];
	std::cout << "BOTTOM LEFT POINT: " << order.p0 << "\n";
	//pop bottom left point, sort, put bottom left point in front
	points.erase(std::remove(points.begin(), points.end(), order.p0), points.end());
	std::sort(points.begin(), points.end(), order);
	points.insert(points.begin(), order.p0);
	std::cout << "LENGTH SORTED POINTS: " << points.size() << "\n";
	std::cout << "FIRST 5 SORTED POINTS: " << "\n";
	std::cout << points[0] << "\n";
	std::cout << points[1] << "\n";
	std::cout << points[2] << "\n";
	std::cout << points[3] << "\n";
	std::cout << points[4] << "\n";
	Polygon hull;
	// TODO
	//iterate through sorted points (and prev/next) to determine which groups of three 
	//form right turns, if so remove middle (current) point
	hull.insert(hull.end(), points[0]); 
	hull.insert(hull.end(), points[1]); 
	hull.insert(hull.end(), points[2]); 

	for (int i = 3; i < points.size(); i++){
		while (hull.size() > 1 && salientAngle(hull[hull.size()-2], hull[hull.size()-1], points[i]) == true){
			hull.pop_back();
		}
		hull.push_back(points[i]);
	}

	return hull;
}

////////////////////////////////////////////////////////////////////////////////

//task 4
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

	return points;
}

void save_obj(const std::string &filename, Polygon &poly) {
	std::ofstream out(filename);
	if (!out.is_open()) {
		throw std::runtime_error("failed to open file " + filename);
	}
	out << std::fixed;
	for (const auto &v : poly) {
		out << "v " << v.real() << ' ' << v.imag() << " 0\n";
	}
	for (size_t i = 0; i < poly.size(); ++i) {
		out << "l " << i+1 << ' ' << 1+(i+1)%poly.size() << "\n";
	}
	out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char * argv[]) {
	if (argc <= 2) {
		std::cerr << "Usage: " << argv[0] << " points.xyz output.obj" << std::endl;
	}
	std::vector<Point> points = load_xyz("/Users/erinmcgowan/Archive/points.xyz");
	Polygon hull = convex_hull(points);
	std::cout << "LENGTH OF CONVEX HULL VECTOR: " << hull.size() << "\n";
	std::cout << "HULL SAVED IN: "<< argv[2] << "\n";
	save_obj(argv[2], hull);
	return 0;
}
