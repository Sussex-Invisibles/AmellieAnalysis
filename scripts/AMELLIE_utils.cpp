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
triangle::triangle(double x_a, double x_b, double x_c, double y_a, double y_b, double y_c) {
    points[0] = x_a; points[1] = x_b; points[2] = x_c;
    points[3] = y_a; points[4] = y_b; points[5] = y_c;
    // Twice the area of the triangle (twice to save on computation later)
    Area2 = abs(points[0] * (points[4] - points[5]) + points[1] * (points[5] - points[3])
                + points[2] * (points[3] - points[4]));
}

// member functions

// Get/Set triangle vertex coordinates
double& triangle::X_a() {return points[0];}
double& triangle::X_b() {return points[1];}
double& triangle::X_c() {return points[2];}
double& triangle::Y_a() {return points[3];}
double& triangle::Y_b() {return points[4];}
double& triangle::Y_c() {return points[5];}

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
    if (Area2 == 0.0) {
        std::cout << "Triangle area is zero." << std::endl;
        if (point_x == points[0] and point_y == points[3]) {
            return lim_count;  // still might be on the point of the triangle
        } else {
            return false;
        }
    }

    // Calculate (twice) the area of each triangle formed by the point, and check if it is zero successively
    double A_a = abs(point_x * (points[4] - points[5]) + points[1] * (points[5] - point_y)
                    + points[2] * (point_y - points[4]));
    if (A_a == 0.0) {return lim_count;} // if area is zero, point is on edge of triangle.

    double A_b = abs(points[0] * (point_y - points[5]) + point_x * (points[5] - points[3])
                    + points[2] * (points[3] - point_y));
    if (A_b == 0.0) {return lim_count;} // if area is zero, point is on edge of triangle.

    double A_c = abs(points[0] * (points[4] - point_y) + points[1] * (point_y - points[3])
                    + point_x * (points[3] - points[4]));
    if (A_c == 0.0) {return lim_count;} // if area is zero, point is on edge of triangle.

    // Check sum of areas of point triangles adds up to original triangle aread (tolerance of 0.1%).
    // If they equal, the point is inside triangle. If the sum is larger, the point is outside.
    double frac_diff = ((A_a + A_b + A_c) / Area2) - 1.0;
    if (frac_diff <= 0.001) {
        return true;
    } else {
        return false;
    }
}


/* ~~~~~~~~~ operator overload ~~~~~~~~~ */
double& triangle::operator [] (int i) {
    assert((i < 6 && i >= 0) && "Index out of range");
    return points[i];
}