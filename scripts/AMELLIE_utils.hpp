// header guard:
#ifndef AMELLIE_utils
#define AMELLIE_utils

// include
#include <iostream>

class triangle {
    private:
        double x_a; double x_b; double x_c;
        double y_a; double y_b; double y_c;

    public:
        // constructor
        triangle(double Xa, double Xb, double Xc, double Ya, double Yb, double Yc);

        // member functions
        double& X_a();
        double& X_b();
        double& X_c();
        double& Y_a();
        double& Y_b();
        double& Y_c();

        // double X_a(double d);
        // double X_b(double d);
        // double X_c(double d);
        // double Y_a(double d);
        // double Y_b(double d);
        // double Y_c(double d);

        bool check_point_inside_triangle(const double point_x, const double point_y, const bool lim_count=true);
};

//end header guard
#endif