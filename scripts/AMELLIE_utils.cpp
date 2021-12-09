#include "AMELLIE_utils.hpp"


/**
 * @brief Construct a new triangle::triangle object
 * 
 * @param x_a 
 * @param x_b 
 * @param x_c 
 * @param y_a 
 * @param y_b 
 * @param y_c 
 */
triangle::triangle(double Xa, double Xb, double Xc, double Ya, double Yb, double Yc) {
    x_a = Xa; x_b = Xb; x_c = Xc; y_a = Ya; y_b = Yb; y_c = Yc;
}

// member functions

// Get/Set triangle vertex coordinates
double& triangle::X_a() {return x_a;}
double& triangle::X_b() {return x_b;}
double& triangle::X_c() {return x_c;}
double& triangle::Y_a() {return y_a;}
double& triangle::Y_b() {return y_b;}
double& triangle::Y_c() {return y_c;}

// // Set triangle vertex coordinates
// double triangle::X_a(double d) {x_a = d;}
// double triangle::X_b(double d) {x_b = d;}
// double triangle::X_c(double d) {x_c = d;}
// double triangle::Y_a(double d) {y_a = d;}
// double triangle::Y_b(double d) {y_b = d;}
// double triangle::Y_c(double d) {y_c = d;}

/**
 * @brief Checks if point is inside triangle by using the sum of triangles method:
 * if the sum of the area of the three triangles formed by the point and two vertices of the triangle add up to
 * the area of the original triangle, the point is inside. It is outside if the sum is larger. 
 * 
 * @param point_x 
 * @param point_y 
 * @param lim_count true (default) = point on triangle edge counts as inside, false = point on edge counts as outside.
 * @return true 
 * @return false 
 */
bool triangle::check_point_inside_triangle(const double point_x, const double point_y, const bool lim_count) {
    // Calculate (twice, since the factor 0.5 is irrelevent here) the triangle's area
    double A = abs(x_a * (y_b - y_c) + x_b * (y_c - y_a) + x_c * (y_a - y_b));

    if (A == 0.0) {
        std::cout << "Triangle area is zero." << std::endl;
        if (point_x == x_a and point_y == y_a) {
            return lim_count;  // still might be on the point of the triangle
        } else {
            return false;
        }
    }

    // Calculate (twice) the area of each triangle formed by the point, and check if it is zero successively
    double A_a = abs(point_x * (y_b - y_c) + x_b * (y_c - point_y) + x_c * (point_y - y_b));
    if (A_a == 0.0) {return lim_count;} // if area is zero, point is on edge of triangle.

    double A_b = abs(x_a * (point_y - y_c) + point_x * (y_c - y_a) + x_c * (y_a - point_y));
    if (A_b == 0.0) {return lim_count;} // if area is zero, point is on edge of triangle.

    double A_c = abs(x_a * (y_b - point_y) + x_b * (point_y - y_a) + point_x * (y_a - y_b));
    if (A_c == 0.0) {return lim_count;} // if area is zero, point is on edge of triangle.

    // Check sum of areas of point triangles adds up to original triangle aread (tolerance of 0.1%).
    // If they equal, the point is inside triangle. If the sum is larger, the point is outside.
    double frac_diff = ((A_a + A_b + A_c) / A) - 1.0;
    if (frac_diff <= 0.001) {
        return true;
    } else {
        return false;
    }
}