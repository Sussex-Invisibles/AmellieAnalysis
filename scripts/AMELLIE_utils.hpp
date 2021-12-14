// header guard:
#ifndef AMELLIE_utils
#define AMELLIE_utils

// include
#include <iostream>
#include <fstream>
#include <assert.h>
#include <algorithm>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TObject.h>
#include <TList.h>
#include <TVectorD.h>

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


class HistList {
    private:
        unsigned int length;
        std::vector<TH2F*> tracking_hists;
        std::vector<TH2F*> region_hists;
        std::vector<TH2F*> direct_hists;
        std::vector<TH2F*> reflected_hists;

    public:
    // constructor
    HistList(std::string tracking_file, std::vector<std::string> name_list
                = {"hReemissionResTimeVsCosTheta", "hPmtResTimeVsCosTheta", "hNoiseResTimeVsCosTheta", 
                "hSingleScatterResTimeVsCosTheta", "hOtherEffectResTimeVsCosTheta", "hNoEffectResTimeVsCosTheta", 
                "hNearReflectResTimeVsCosTheta", "hRopesResTimeVsCosTheta", "hPMTReflectionResTimeVsCosTheta", 
                "hExtWaterScatterResTimeVsCosTheta", "hInnerAvReflectionResTimeVsCosTheta", "hMultipleEffectResTimeVsCosTheta",
                "hAVPipesResTimeVsCosTheta", "hAcrylicScatterResTimeVsCosTheta", "OtherScatterResTimeVsCosTheta"});

    // return suitable histograms
    unsigned int len();
    std::vector<TH2F*>& Tracking_Hists();
    std::vector<TH2F*>& Region_Hists();
    std::vector<TH2F*>& Direct_Hists();
    std::vector<TH2F*>& Reflected_Hists();
};

//end header guard
#endif