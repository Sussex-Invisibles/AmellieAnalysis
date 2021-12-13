// header guard:
#ifndef AMELLIE_utils
#define AMELLIE_utils

// include
#include <iostream>
#include <assert.h>

class triangle {
    private:
        double points[6]; //x_a, x_b, x_c, y_a, y_b, y_c
        double Area2;

    public:
        // constructors
        triangle(double x_a, double x_b, double x_c, double y_a, double y_b, double y_c);

        // member functions
        double& X_a();
        double& X_b();
        double& X_c();
        double& Y_a();
        double& Y_b();
        double& Y_c();

        bool check_point_inside_triangle(const double point_x, const double point_y, const bool lim_count=true);

    /* ~~~~~~~~~ operator overload ~~~~~~~~~ */

        // access elements easier
        double& operator [] (int i);
};

//end header guard
#endif