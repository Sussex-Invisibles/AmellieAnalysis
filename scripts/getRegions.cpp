//Compile: g++ -g -std=c++1y -o getRegions.exe getRegions.cpp `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux
/*
Rewrite of the region selection code, after conversations with Lisa.
Create plots of the phase space showing the relationship between different parameters [x_a, y_a, x_b_c, y_b, y_c]
*/
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <TCanvas.h>
#include <TLine.h>
#include "AMELLIE_utils.hpp"

int CalculateRegions(std::string inputFile, int nbins);
int OptimiseDivideAndConquer(std::string inputFile, int nbins, bool verbose, bool debug, bool extraInfo, std::string signal_param);
std::vector<double> GetThreePoints(double bestPoint, double worstPoint, std::vector<double> originalPoints);
std::vector<double> GetFOMs(std::vector<double> points, std::vector<double> fixedPoints, int numVar, TH2F *allPathsHist, TH2F *reEmittedHist, TH2F *scatteredHist, std::string signal);
std::vector<double> GetBestFOM(std::vector<double> FOMs, std::vector<double> points);
std::vector<TH2F*> GetRegionSelectedHists(std::vector<double> finalPoints, TH2F *hReEmittedPaths, TH2F *hAllPaths, TH2F *hNoisePaths, TH2F *hSingleScatterPaths, TH2F *hOtherPaths, TH2F *hNoEffectPaths, TH2F *hNearReflectPaths, TH2F *hRopesPaths, TH2F *hPMTReflectionPaths, TH2F *hExtWaterScatterPaths, TH2F *hInnerAvReflectPaths, TH2F *hMultipleEffectPaths, TH2F *hAVPipesPaths, TH2F *hAcrylicPaths, TH2F *hOtherScatterPaths);
std::vector<double> CheckPoints(std::vector<double> points, std::vector<double> fixedPoints, int numVar);

int main(int argc, char** argv){
    std::string file = argv[1];
    int nbins = std::stoi(argv[2]);
    bool verbose = std::stoi(argv[3]);
    bool debug = std::stoi(argv[4]);
    bool extraInfo = std::stoi(argv[5]);
    std::string signal = argv[6];
    auto t1 = std::chrono::high_resolution_clock::now();
    int status = OptimiseDivideAndConquer(file, nbins, verbose, debug, extraInfo, signal); 
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "Script took " << std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() / 1E6 << " s to execute" << std::endl;

    return 0;
}

int OptimiseDivideAndConquer(std::string inputFile, int nbins, bool verbose, bool debug, bool extraInfo, std::string signal_param){

    auto timeStart = std::chrono::high_resolution_clock::now();

    if(signal_param != "reemitted" and signal_param != "scattered" and signal_param != "attenuated"){ // attenuated = scattered + reemitted
        std::cout << "Invalid signal! Use: reemitted, scattered, attenuated";
        return 1;
    }

    //create histograms and graphs

    TH2F *hReEmittedPaths;
    TH2F *hAllPaths;
    TH2F *hNoisePaths;
    TH2F *hSingleScatterPaths;
    TH2F *hOtherPaths;
    TH2F *hNoEffectPaths;
    TH2F *hNearReflectPaths;
    TH2F *hRopesPaths;
    TH2F *hPMTReflectionPaths;
    TH2F *hExtWaterScatterPaths;
    TH2F *hInnerAvReflectPaths;
    TH2F *hMultipleEffectPaths;
    TH2F *hAVPipesPaths;
    TH2F *hAcrylicPaths;
    TH2F *hOtherScatterPaths;
    TGraph *hFOMx_a = new TGraph();
    TGraph *hFOMx_b = new TGraph();
    TGraph *hFOMx_c = new TGraph();
    TGraph *hFOMy_a = new TGraph();
    TGraph *hFOMy_b = new TGraph();
    TGraph *hFOMy_c = new TGraph();
    TGraph *hFOMmainloop = new TGraph();
    TGraph *hDiffx_a = new TGraph();
    TGraph *hDiffx_b = new TGraph();
    TGraph *hDiffx_c = new TGraph();
    TGraph *hDiffy_a = new TGraph();
    TGraph *hDiffy_b = new TGraph();
    TGraph *hDiffy_c = new TGraph();
    TGraph *hPointx_a = new TGraph();
    TGraph *hPointx_b = new TGraph();
    TGraph *hPointx_c = new TGraph();
    TGraph *hPointy_a = new TGraph();
    TGraph *hPointy_b = new TGraph();
    TGraph *hPointy_c = new TGraph();

    //read in hist file

    TFile *raw_file;

    try{
        raw_file = TFile::Open(inputFile.c_str());
    }
    catch(...){
        std::cout << "Could not open input file." << std::endl;
        return 1;
    }

    //assign hists

    raw_file->GetObject("hReemissionResTimeVsCosTheta",hReEmittedPaths);
    raw_file->GetObject("hPmtResTimeVsCosTheta",hAllPaths);
    raw_file->GetObject("hNoiseResTimeVsCosTheta",hNoisePaths);
    raw_file->GetObject("hSingleScatterResTimeVsCosTheta",hSingleScatterPaths);
    raw_file->GetObject("hOtherEffectResTimeVsCosTheta",hOtherPaths);
    raw_file->GetObject("hNoEffectResTimeVsCosTheta",hNoEffectPaths);
    raw_file->GetObject("hNearReflectResTimeVsCosTheta",hNearReflectPaths);
    raw_file->GetObject("hRopesResTimeVsCosTheta",hRopesPaths);
    raw_file->GetObject("hPMTReflectionResTimeVsCosTheta",hPMTReflectionPaths);
    raw_file->GetObject("hExtWaterScatterResTimeVsCosTheta",hExtWaterScatterPaths);
    raw_file->GetObject("hInnerAvReflectionResTimeVsCosTheta",hInnerAvReflectPaths);
    raw_file->GetObject("hMultipleEffectResTimeVsCosTheta",hMultipleEffectPaths);
    raw_file->GetObject("hAVPipesResTimeVsCosTheta",hAVPipesPaths);
    raw_file->GetObject("hAcrylicScatterResTimeVsCosTheta",hAcrylicPaths);
    raw_file->GetObject("OtherScatterResTimeVsCosTheta",hOtherScatterPaths);

    //set up constants

    double countCalculations = 0;
    double countCycles = 0;
    double x_a_min = -1;
    double x_a_max = 1;
    double x_a_tolerance = (x_a_max - x_a_min) / nbins;
    double x_b_min = -1;
    double x_b_max = 1;
    double x_b_tolerance = (x_b_max - x_b_min) / nbins;
    double x_c_min = -1;
    double x_c_max = 1;
    double x_c_tolerance = (x_c_max - x_c_min) / nbins;
    double y_a_min = hReEmittedPaths->GetYaxis()->GetXmin();
    double y_a_max = hReEmittedPaths->GetYaxis()->GetXmax();
    double y_a_tolerance = (y_a_max - y_a_min) / nbins;
    double y_b_min = hReEmittedPaths->GetYaxis()->GetXmin();
    double y_b_max = hReEmittedPaths->GetYaxis()->GetXmax();
    double y_b_tolerance = (y_b_max - y_b_min) / nbins;
    double y_c_min = hReEmittedPaths->GetYaxis()->GetXmin();
    double y_c_max = hReEmittedPaths->GetYaxis()->GetXmax();
    double y_c_tolerance = (y_c_max - y_c_min) / nbins;

    int nBinsX = hReEmittedPaths->GetXaxis()->GetNbins();
    int nBinsY = hReEmittedPaths->GetYaxis()->GetNbins();

    double xBinWidth = hReEmittedPaths->GetXaxis()->GetBinCenter(2) - hReEmittedPaths->GetXaxis()->GetBinCenter(1);
    double yBinWidth = hReEmittedPaths->GetYaxis()->GetBinCenter(2) - hReEmittedPaths->GetYaxis()->GetBinCenter(1);

    if(x_a_tolerance < xBinWidth){
        x_a_tolerance = xBinWidth;
        std::cout << "x_a_tolerance is set to bin width, nbins: " << (x_a_max - x_a_min) / xBinWidth << std::endl;
    }
    if(x_b_tolerance < xBinWidth){
        x_b_tolerance = xBinWidth;
        std::cout << "x_b_tolerance is set to bin width, nbins: " << (x_b_max - x_b_min) / xBinWidth << std::endl;
    }
    if(x_c_tolerance < xBinWidth){
        x_c_tolerance = xBinWidth;
        std::cout << "x_c_tolerance is set to bin width, nbins: " << (x_c_max - x_c_min) / xBinWidth << std::endl;
    }
    if(y_a_tolerance < yBinWidth){
        y_a_tolerance = yBinWidth;
        std::cout << "y_a_tolerance is set to bin width, nbins: " << (y_a_max - y_a_min) / yBinWidth << std::endl;
    }
    if(y_b_tolerance < yBinWidth){
        y_b_tolerance = yBinWidth;
        std::cout << "y_b_tolerance is set to bin width, nbins: " << (y_b_max - y_b_min) / yBinWidth << std::endl;
    }
    if(y_c_tolerance < yBinWidth){
        y_c_tolerance = yBinWidth;
        std::cout << "y_c_tolerance is set to bin width, nbins: " << (y_c_max - y_c_min) / yBinWidth << std::endl;
    }

    if(debug) std::cout << "x_a_tolerance: " << x_a_tolerance << ", x_b_tolerance: " << x_b_tolerance << std::endl;

    double x_a_diff = 9999999999;
    double x_b_diff = 9999999999;
    double x_c_diff = 9999999999;
    double y_a_diff = 9999999999;
    double y_b_diff = 9999999999;
    double y_c_diff = 9999999999;

    double x_a_temp_diff = 9999999999;
    double x_b_temp_diff = 9999999999;
    double x_c_temp_diff = 9999999999;
    double y_a_temp_diff = 9999999999;
    double y_b_temp_diff = 9999999999;
    double y_c_temp_diff = 9999999999;

    auto timeSetup = std::chrono::high_resolution_clock::now();

    if(verbose) std::cout << "Time to read in file and set up hists: " << std::chrono::duration_cast<std::chrono::microseconds>(timeSetup - timeStart).count() / 1E6 << "s" << std::endl;
    // do loop: x_a, x_bc, y_a, y_b, y_c
    // points: a vector containing the points to change for three regions
    // fixedPoints: a vector containing all the static points

    auto timeStartLoop = std::chrono::high_resolution_clock::now();

    std::vector<double> fixedPoints = {x_a_min, x_b_max, x_c_max, y_a_min, y_b_max, y_c_min};
    bool firstRun = true;

    int numMainLoopIterations = 0;
    int numx_a_iterations = 0;
    int numx_b_iterations = 0;
    int numx_c_iterations = 0;
    int numy_a_iterations = 0;
    int numy_b_iterations = 0;
    int numy_c_iterations = 0;

    double x_a_prev_best_point_main;
    double x_b_prev_best_point_main;
    double x_c_prev_best_point_main;
    double y_a_prev_best_point_main;
    double y_b_prev_best_point_main;
    double y_c_prev_best_point_main;

    while(((x_a_diff > x_a_tolerance or x_b_diff > x_b_tolerance or x_c_diff > x_c_tolerance or y_a_diff > y_a_tolerance or y_b_diff > y_b_tolerance or y_c_diff > y_c_tolerance) and !firstRun) or firstRun){
        bool first_x_a_run = true;
        bool first_x_b_run = true;
        bool first_x_c_run = true;
        bool first_y_a_run = true;
        bool first_y_b_run = true;
        bool first_y_c_run = true;

        double prevBestFOMPoint_x_a = 9999999999;
        double prevBestFOMPoint_x_b = 9999999999;
        double prevBestFOMPoint_x_c = 9999999999;
        double prevBestFOMPoint_y_a = 9999999999;
        double prevBestFOMPoint_y_b = 9999999999;
        double prevBestFOMPoint_y_c = 9999999999;

        std::vector<double> points_x_a;
        std::vector<double> points_x_b;
        std::vector<double> points_x_c;
        std::vector<double> points_y_a;
        std::vector<double> points_y_b;
        std::vector<double> points_y_c;
        std::vector<double> bestworstFOMPoints_x_a;
        std::vector<double> bestworstFOMPoints_x_b;
        std::vector<double> bestworstFOMPoints_x_c;
        std::vector<double> bestworstFOMPoints_y_a;
        std::vector<double> bestworstFOMPoints_y_b;
        std::vector<double> bestworstFOMPoints_y_c;

        double x_a_temp_diff = 9999999999;
        double x_b_temp_diff = 9999999999;
        double x_c_temp_diff = 9999999999;
        double y_a_temp_diff = 9999999999;
        double y_b_temp_diff = 9999999999;
        double y_c_temp_diff = 9999999999;

        if(debug) std::cout << "In main while loop with x_a_temp_diff: " << std::abs(x_a_temp_diff) << std::endl;

        while(std::abs(x_a_temp_diff) > x_a_tolerance){
            if(first_x_a_run){
                points_x_a = {x_a_min, (x_a_max + x_a_min)/2., x_a_max}; //l, c, r
                points_x_a = CheckPoints(points_x_a, fixedPoints, 0);
                first_x_a_run = false;
                if(debug) std::cout << "Set up first x_a_run points: " << points_x_a.at(0) << ", " << points_x_a.at(1) << ", " << points_x_a.at(2) << std::endl;
            }
            std::vector<double> FOMs_x_a = GetFOMs(points_x_a, fixedPoints, 0, hAllPaths, hReEmittedPaths, hSingleScatterPaths, signal_param);
            if(debug) std::cout << "Got FOMs " << FOMs_x_a.at(0) << ", " << FOMs_x_a.at(1) << ", " << FOMs_x_a.at(2) << std::endl;
            
            bestworstFOMPoints_x_a = GetBestFOM(FOMs_x_a, points_x_a);
            if(debug) std::cout << "Got best FOMs" << std::endl;
            
            points_x_a = GetThreePoints(bestworstFOMPoints_x_a.at(0), bestworstFOMPoints_x_a.at(1), points_x_a);
            points_x_a = CheckPoints(points_x_a, fixedPoints, 0);
            if(debug) std::cout << "Got new points: " << points_x_a.at(0) << ", " << points_x_a.at(1) << ", " << points_x_a.at(2) << std::endl;
            
            x_a_temp_diff = std::abs(points_x_a.at(1) - prevBestFOMPoint_x_a);
            if(debug) std::cout << "Got difference: " << x_a_temp_diff << std::endl;
            
            prevBestFOMPoint_x_a = points_x_a.at(1);
            if(extraInfo) hFOMx_a->SetPoint(numx_a_iterations, numx_a_iterations, bestworstFOMPoints_x_a.at(2));
            numx_a_iterations++;
        }

        if(debug) std::cout << "Setting fixed point: " << bestworstFOMPoints_x_a.at(0) << std::endl;
        fixedPoints.at(0) = bestworstFOMPoints_x_a.at(0);
        if(debug) std::cout << "" << std::endl;

        while(std::abs(x_b_temp_diff) > x_b_tolerance){
            if(first_x_b_run){
                points_x_b = {x_b_min, (x_b_max + x_b_min)/2., x_b_max}; //l, c, r
                points_x_b = CheckPoints(points_x_b, fixedPoints, 1);
                first_x_b_run = false;
                if(debug) std::cout << "First x_bc_run points: " << points_x_b.at(0) << ", " << points_x_b.at(1) << ", " << points_x_b.at(2) << std::endl;
            }
            std::vector<double> FOMs_x_b = GetFOMs(points_x_b, fixedPoints, 1, hAllPaths, hReEmittedPaths, hSingleScatterPaths, signal_param);
            if(debug) std::cout << "Got FOMs: " << FOMs_x_b.at(0) << ", " << FOMs_x_b.at(1) << ", " << FOMs_x_b.at(2) << std::endl;
            
            bestworstFOMPoints_x_b = GetBestFOM(FOMs_x_b, points_x_b);
            if(debug) std::cout << "Got best FOMs" << std::endl;
            
            points_x_b = GetThreePoints(bestworstFOMPoints_x_b.at(0), bestworstFOMPoints_x_b.at(1), points_x_b);
            points_x_b = CheckPoints(points_x_b, fixedPoints, 1);
            if(debug) std::cout << "Got new points: " << points_x_b.at(0) << ", " << points_x_b.at(1) << ", " << points_x_b.at(2) << std::endl;
            
            x_b_temp_diff = std::abs(points_x_b.at(1) - prevBestFOMPoint_x_b);
            if(debug) std::cout << "Got difference: " << x_b_temp_diff << std::endl;
            
            prevBestFOMPoint_x_b = points_x_b.at(1);
            if(extraInfo) hFOMx_b->SetPoint(numx_b_iterations, numx_b_iterations, bestworstFOMPoints_x_b.at(2));
            numx_b_iterations++;
        }
        if(debug) std::cout << "Setting fixed point: " << bestworstFOMPoints_x_b.at(0) << std::endl;
        fixedPoints.at(1) = bestworstFOMPoints_x_b.at(0);
        if(debug) std::cout << "" << std::endl;

        while(std::abs(x_c_temp_diff) > x_c_tolerance){
            if(first_x_c_run){
                points_x_c = {x_c_min, (x_c_max + x_c_min)/2., x_c_max}; //l, c, r
                points_x_c = CheckPoints(points_x_c, fixedPoints, 2);
                first_x_c_run = false;
                if(debug) std::cout << "First x_c_run points: " << points_x_c.at(0) << ", " << points_x_c.at(1) << ", " << points_x_c.at(2) << std::endl;
            }
            std::vector<double> FOMs_x_c = GetFOMs(points_x_c, fixedPoints, 2, hAllPaths, hReEmittedPaths, hSingleScatterPaths, signal_param);
            if(debug) std::cout << "Got FOMs: " << FOMs_x_c.at(0) << ", " << FOMs_x_c.at(1) << ", " << FOMs_x_c.at(2) << std::endl;
            
            bestworstFOMPoints_x_c = GetBestFOM(FOMs_x_c, points_x_c);
            if(debug) std::cout << "Got best FOMs" << std::endl;
            
            points_x_c = GetThreePoints(bestworstFOMPoints_x_c.at(0), bestworstFOMPoints_x_c.at(1), points_x_c);
            points_x_c = CheckPoints(points_x_c, fixedPoints, 2);
            if(debug) std::cout << "Got new points: " << points_x_c.at(0) << ", " << points_x_c.at(1) << ", " << points_x_c.at(2) << std::endl;
            
            x_c_temp_diff = std::abs(points_x_c.at(1) - prevBestFOMPoint_x_c);
            if(debug) std::cout << "Got difference: " << x_c_temp_diff << std::endl;
            
            prevBestFOMPoint_x_c = points_x_c.at(1);
            if(extraInfo) hFOMx_c->SetPoint(numx_c_iterations, numx_c_iterations, bestworstFOMPoints_x_c.at(2));
            numx_c_iterations++;
        }
        if(debug) std::cout << "Setting fixed point: " << bestworstFOMPoints_x_c.at(0) << std::endl;
        fixedPoints.at(2) = bestworstFOMPoints_x_c.at(0);
        if(debug) std::cout << "" << std::endl;

        while(std::abs(y_a_temp_diff) > y_a_tolerance){
            if(first_y_a_run){
                if(debug) std::cout << "Setting up first y_a_run" << std::endl;
                points_y_a = {y_a_min, (y_a_max + y_a_min)/2., y_a_max}; //l, c, r
                points_y_a = CheckPoints(points_y_a, fixedPoints, 3);
                first_y_a_run = false;
                if(debug) std::cout << "Set up first y_a_run points: " << points_y_a.at(0) << ", " << points_y_a.at(1) << ", " << points_y_a.at(2) << std::endl;
            }
            std::vector<double> FOMs_y_a = GetFOMs(points_y_a, fixedPoints, 3, hAllPaths, hReEmittedPaths, hSingleScatterPaths, signal_param);
            if(debug) std::cout << "Got FOMs " << FOMs_y_a.at(0) << ", " << FOMs_y_a.at(1) << ", " << FOMs_y_a.at(2) << std::endl;
            
            bestworstFOMPoints_y_a = GetBestFOM(FOMs_y_a, points_y_a);
            if(debug) std::cout << "Got best FOMs" << std::endl;
            
            points_y_a = GetThreePoints(bestworstFOMPoints_y_a.at(0), bestworstFOMPoints_y_a.at(1), points_y_a);
            points_y_a = CheckPoints(points_y_a, fixedPoints, 3);
            if(debug) std::cout << "Got new points: " << points_y_a.at(0) << ", " << points_y_a.at(1) << ", " << points_y_a.at(2) << std::endl;
            
            y_a_temp_diff = std::abs(points_y_a.at(1) - prevBestFOMPoint_y_a);
            if(debug) std::cout << "Got difference: " << y_a_temp_diff << std::endl;
            
            prevBestFOMPoint_y_a = points_y_a.at(1);
            if(extraInfo) hFOMy_a->SetPoint(numy_a_iterations, numy_a_iterations, bestworstFOMPoints_y_a.at(2));
            numy_a_iterations++;
        }
        if(debug) std::cout << "Setting fixed point: " << bestworstFOMPoints_y_a.at(0) << std::endl;
        fixedPoints.at(3) = bestworstFOMPoints_y_a.at(0);
        if(debug) std::cout << "" << std::endl;

        while(std::abs(y_b_temp_diff) > y_b_tolerance){
            if(first_y_b_run){
                if(debug) std::cout << "Setting up first y_b_run" << std::endl;
                points_y_b = {y_b_min, (y_b_max + y_b_min)/2., y_b_max}; //l, c, r
                points_y_b = CheckPoints(points_y_b, fixedPoints, 4);
                first_y_b_run = false;
                if(debug) std::cout << "Set up first y_b_run points: " << points_y_b.at(0) << ", " << points_y_b.at(1) << ", " << points_y_b.at(2) << std::endl;
            }
            std::vector<double> FOMs_y_b = GetFOMs(points_y_b, fixedPoints, 4, hAllPaths, hReEmittedPaths, hSingleScatterPaths, signal_param);
            if(debug) std::cout << "Got FOMs " << FOMs_y_b.at(0) << ", " << FOMs_y_b.at(1) << ", " << FOMs_y_b.at(2) << std::endl;
            
            bestworstFOMPoints_y_b = GetBestFOM(FOMs_y_b, points_y_b);
            if(debug) std::cout << "Got best FOMs" << std::endl;
            
            points_y_b = GetThreePoints(bestworstFOMPoints_y_b.at(0), bestworstFOMPoints_y_b.at(1), points_y_b);
            points_y_b = CheckPoints(points_y_b, fixedPoints, 4);
            if(debug) std::cout << "Got new points: " << points_y_b.at(0) << ", " << points_y_b.at(1) << ", " << points_y_b.at(2) << std::endl;
            
            y_b_temp_diff = std::abs(points_y_b.at(1) - prevBestFOMPoint_y_b);
            if(debug) std::cout << "Got difference: " << y_b_temp_diff << std::endl;
            
            prevBestFOMPoint_y_b = points_y_b.at(1);
            if(extraInfo) hFOMy_b->SetPoint(numy_b_iterations, numy_b_iterations, bestworstFOMPoints_y_b.at(2));
            numy_b_iterations++;
        }
        if(debug) std::cout << "Setting fixed point: " << bestworstFOMPoints_y_b.at(0) << std::endl;
        fixedPoints.at(4) = bestworstFOMPoints_y_b.at(0);
        if(debug) std::cout << "" << std::endl;

        while(std::abs(y_c_temp_diff) > y_c_tolerance){
            if(first_y_c_run){
                if(debug) std::cout << "Setting up first y_c_run" << std::endl;
                points_y_c = {y_c_min, (y_c_max + y_c_min)/2., y_c_max}; //l, c, r
                points_y_c = CheckPoints(points_y_c, fixedPoints, 5);
                first_y_c_run = false;
                if(debug) std::cout << "Set up first y_c_run points: " << points_y_c.at(0) << ", " << points_y_c.at(1) << ", " << points_y_c.at(2) << std::endl;
            }
            std::vector<double> FOMs_y_c = GetFOMs(points_y_c, fixedPoints, 5, hAllPaths, hReEmittedPaths, hSingleScatterPaths, signal_param);
            if(debug) std::cout << "Got FOMs " << FOMs_y_c.at(0) << ", " << FOMs_y_c.at(1) << ", " << FOMs_y_c.at(2) << std::endl;
            
            bestworstFOMPoints_y_c = GetBestFOM(FOMs_y_c, points_y_c);
            if(debug) std::cout << "Got best FOMs" << std::endl;
            
            points_y_c = GetThreePoints(bestworstFOMPoints_y_c.at(0), bestworstFOMPoints_y_c.at(1), points_y_c);
            points_y_c = CheckPoints(points_y_c, fixedPoints, 5);
            if(debug) std::cout << "Got new points: " << points_y_c.at(0) << ", " << points_y_c.at(1) << ", " << points_y_c.at(2) << std::endl;
            
            y_c_temp_diff = std::abs(points_y_c.at(1) - prevBestFOMPoint_y_c);
            if(debug) std::cout << "Got difference: " << y_c_temp_diff << std::endl;
            
            prevBestFOMPoint_y_c = points_y_c.at(1);
            if(extraInfo) hFOMy_c->SetPoint(numy_c_iterations, numy_c_iterations, bestworstFOMPoints_y_c.at(2));
            numy_c_iterations++;
        }
        if(debug) std::cout << "Setting fixed point: " << bestworstFOMPoints_y_c.at(0) << std::endl;
        fixedPoints.at(5) = bestworstFOMPoints_y_c.at(0);
        if(debug) std::cout << "" << std::endl;

        if(debug) std::cout << "Getting main differences" << std::endl;

        if(x_a_diff == 9999999999){
            x_a_diff = std::abs(bestworstFOMPoints_x_a.at(0));
            x_a_prev_best_point_main = bestworstFOMPoints_x_a.at(0);
        }
        else{
            x_a_diff = std::abs(x_a_prev_best_point_main - bestworstFOMPoints_x_a.at(0));
            x_a_prev_best_point_main = bestworstFOMPoints_x_a.at(0);
            firstRun = false;
        }

        if(x_b_diff == 9999999999){
            x_b_diff = std::abs(bestworstFOMPoints_x_b.at(0));
            x_b_prev_best_point_main = bestworstFOMPoints_x_b.at(0);
        }
        else{
            x_b_diff = std::abs(x_b_prev_best_point_main - bestworstFOMPoints_x_b.at(0));
            x_b_prev_best_point_main = bestworstFOMPoints_x_b.at(0);
        }

        if(x_c_diff == 9999999999){
            x_c_diff = std::abs(bestworstFOMPoints_x_c.at(0));
            x_c_prev_best_point_main = bestworstFOMPoints_x_c.at(0);
        }
        else{
            x_c_diff = std::abs(x_c_prev_best_point_main - bestworstFOMPoints_x_c.at(0));
            x_c_prev_best_point_main = bestworstFOMPoints_x_c.at(0);
        }

        if(y_a_diff == 9999999999){
            y_a_diff = std::abs(bestworstFOMPoints_y_a.at(0));
            y_a_prev_best_point_main = bestworstFOMPoints_y_a.at(0);
        }
        else{
            y_a_diff = std::abs(y_a_prev_best_point_main - bestworstFOMPoints_y_a.at(0));
            y_a_prev_best_point_main = bestworstFOMPoints_y_a.at(0);
        }

        if(y_b_diff == 9999999999){
            y_b_diff = std::abs(bestworstFOMPoints_y_b.at(0));
            y_b_prev_best_point_main = bestworstFOMPoints_y_b.at(0);
        }
        else{
            y_b_diff = std::abs(y_b_prev_best_point_main - bestworstFOMPoints_y_b.at(0));
            y_b_prev_best_point_main = bestworstFOMPoints_y_b.at(0);
        }

        if(y_c_diff == 9999999999){
            y_c_diff = std::abs(bestworstFOMPoints_y_c.at(0));
            y_c_prev_best_point_main = bestworstFOMPoints_y_c.at(0);
        }
        else{
            y_c_diff = std::abs(y_c_prev_best_point_main - bestworstFOMPoints_y_c.at(0));
            y_c_prev_best_point_main = bestworstFOMPoints_y_c.at(0);
        }
        if(debug) std::cout << "Got main differences" << std::endl;
        if(debug) std::cout << "" << std::endl;
        if(extraInfo) hFOMmainloop->SetPoint(numMainLoopIterations, numMainLoopIterations, bestworstFOMPoints_y_c.at(2));
        if(extraInfo) hDiffx_a->SetPoint(numMainLoopIterations, numMainLoopIterations, x_a_diff);
        if(extraInfo) hDiffx_b->SetPoint(numMainLoopIterations, numMainLoopIterations, x_b_diff);
        if(extraInfo) hDiffx_c->SetPoint(numMainLoopIterations, numMainLoopIterations, x_c_diff);
        if(extraInfo) hDiffy_a->SetPoint(numMainLoopIterations, numMainLoopIterations, y_a_diff);
        if(extraInfo) hDiffy_b->SetPoint(numMainLoopIterations, numMainLoopIterations, y_b_diff);
        if(extraInfo) hDiffy_c->SetPoint(numMainLoopIterations, numMainLoopIterations, y_c_diff);
        if(extraInfo) hPointx_a->SetPoint(numMainLoopIterations, numMainLoopIterations, fixedPoints.at(0));
        if(extraInfo) hPointx_b->SetPoint(numMainLoopIterations, numMainLoopIterations, fixedPoints.at(1));
        if(extraInfo) hPointx_c->SetPoint(numMainLoopIterations, numMainLoopIterations, fixedPoints.at(2));
        if(extraInfo) hPointy_a->SetPoint(numMainLoopIterations, numMainLoopIterations, fixedPoints.at(3));
        if(extraInfo) hPointy_b->SetPoint(numMainLoopIterations, numMainLoopIterations, fixedPoints.at(4));
        if(extraInfo) hPointy_c->SetPoint(numMainLoopIterations, numMainLoopIterations, fixedPoints.at(5));
        numMainLoopIterations++;
    }

    //now use final result to get region, write to file

    if(verbose) std::cout << "Region found" << std::endl;
    if(verbose) std::cout << "    x_a: " << fixedPoints.at(0) << std::endl;
    if(verbose) std::cout << "    x_b: " << fixedPoints.at(1) << std::endl;
    if(verbose) std::cout << "    x_c: " << fixedPoints.at(2) << std::endl;
    if(verbose) std::cout << "    y_a: " << fixedPoints.at(3) << std::endl;
    if(verbose) std::cout << "    y_b: " << fixedPoints.at(4) << std::endl;
    if(verbose) std::cout << "    y_c: " << fixedPoints.at(5) << std::endl;
    if(verbose) std::cout << "There were " << numMainLoopIterations << " iterations of the main loop" << std::endl;

    std::vector<TH2F*> regionSelectedHists = GetRegionSelectedHists(fixedPoints, hReEmittedPaths, hAllPaths, hNoisePaths, hSingleScatterPaths, hOtherPaths, hNoEffectPaths, hNearReflectPaths, hRopesPaths, hPMTReflectionPaths, hExtWaterScatterPaths, hInnerAvReflectPaths, hMultipleEffectPaths, hAVPipesPaths, hAcrylicPaths, hOtherScatterPaths);
    // get file name from path+filename string
    std::size_t botDirPos = inputFile.find_last_of("/");
    std::string filename = inputFile.substr(botDirPos+1, inputFile.length());
    std::string saveroot = "region_selected_hists_" + signal_param + "_" + filename;
    TFile *rootfile = new TFile(saveroot.c_str(),"RECREATE");
    
    rootfile->cd();
    for(int i=0; i<regionSelectedHists.size();i++){
        TH2F *tempHist = regionSelectedHists.at(i);
        tempHist->Write();
    }
    if(extraInfo){
        hFOMx_a->SetName("hFOMx_a");
        hFOMx_b->SetName("hFOMx_b");
        hFOMx_c->SetName("hFOMx_c");
        hFOMy_a->SetName("hFOMy_a");
        hFOMy_b->SetName("hFOMy_b");
        hFOMy_c->SetName("hFOMy_c");
        hFOMmainloop->SetName("hFOMmainloop");
        hDiffx_a->SetName("hDiffx_a");
        hDiffx_b->SetName("hDiffx_b");
        hDiffx_c->SetName("hDiffx_c");
        hDiffy_a->SetName("hDiffy_a");
        hDiffy_b->SetName("hDiffy_b");
        hDiffy_c->SetName("hDiffy_c");
        hPointx_a->SetName("hPointx_a");
        hPointx_b->SetName("hPointx_b");
        hPointx_c->SetName("hPointx_c");
        hPointy_a->SetName("hPointy_a");
        hPointy_b->SetName("hPointy_b");
        hPointy_c->SetName("hPointy_c");
        hFOMx_a->Write();
        hFOMx_b->Write();
        hFOMx_c->Write();
        hFOMy_a->Write();
        hFOMy_b->Write();
        hFOMy_c->Write();
        hFOMmainloop->Write();
        hDiffx_a->Write();
        hDiffx_b->Write();
        hDiffx_c->Write();
        hDiffy_a->Write();
        hDiffy_b->Write();
        hDiffy_c->Write();
        hPointx_a->Write();
        hPointx_b->Write();
        hPointx_c->Write();
        hPointy_a->Write();
        hPointy_b->Write();
        hPointy_c->Write();
    }
    rootfile->Write();
    rootfile->Close();

    return 0;
}

std::vector<double> GetThreePoints(double bestPoint, double worstPoint, std::vector<double> originalPoints){
    std::vector<double> newPoints;
    if(worstPoint == originalPoints.at(0)){ //left point worst
        newPoints = {originalPoints.at(1), (originalPoints.at(1) + originalPoints.at(2)) / 2., originalPoints.at(2)};
    }
    else if(worstPoint == originalPoints.at(2)){ //right point worst
        newPoints = {originalPoints.at(0), (originalPoints.at(0) + originalPoints.at(1)) / 2., originalPoints.at(1)};
    }
    else if(worstPoint == originalPoints.at(1)){ //middle point worst, remove next worse point
        if(originalPoints.at(0) == bestPoint){ //right point worst
            newPoints = {originalPoints.at(0), (originalPoints.at(0) + originalPoints.at(1)) / 2., originalPoints.at(1)};
        }
        else if(originalPoints.at(2) == bestPoint){ //left point worst
            newPoints = {originalPoints.at(1), (originalPoints.at(1) + originalPoints.at(2)) / 2., originalPoints.at(2)};
        }
    }
    else{
        std::cout << "Error! No worst point" <<std::endl;
    }

    return newPoints;
}

std::vector<double> GetFOMs(std::vector<double> points, std::vector<double> fixedPoints, int numVar, TH2F *allPathsHist, TH2F *reEmittedHist, TH2F *scatteredHist, std::string signal){

    //FIXME: Pass in vector of hists?

    double countReEmitted1 = 0;
    double countTotal1 = 0;
    double countReEmitted2 = 0;
    double countTotal2 = 0;
    double countReEmitted3 = 0;
    double countTotal3 = 0;

    //check that values are allowed
    double point1 = points.at(0);
    double point2 = points.at(1);
    double point3 = points.at(2);

    /*
    grad = y2-y1 / x2-x1
    const = y1 - m*x1
    */

    for(int x=1; x<reEmittedHist->GetNbinsX()+1; x++){ //loop over histogram bins
        double xBinCenter = reEmittedHist->GetXaxis()->GetBinCenter(x);
        for(int y=1; y<reEmittedHist->GetNbinsY()+1; y++){
            double yBinCenter = reEmittedHist->GetYaxis()->GetBinCenter(y);
            if(numVar == 0){
                double upperGrad1 = (fixedPoints.at(4) - fixedPoints.at(3)) / (fixedPoints.at(1) - point1);
                double upperIntercept1 = fixedPoints.at(4) - (upperGrad1*fixedPoints.at(1));
                double upperGrad2 = (fixedPoints.at(4) - fixedPoints.at(3)) / (fixedPoints.at(1) - point2);
                double upperIntercept2 = fixedPoints.at(4) - (upperGrad2*fixedPoints.at(1));
                double upperGrad3 = (fixedPoints.at(4) - fixedPoints.at(3)) / (fixedPoints.at(1) - point3);
                double upperIntercept3 = fixedPoints.at(4) - (upperGrad3*fixedPoints.at(1));
                double bottomGrad1 = (fixedPoints.at(5) - fixedPoints.at(3)) / (fixedPoints.at(2) - point1);
                double bottomIntercept1 = fixedPoints.at(5) - (upperGrad1*fixedPoints.at(2));
                double bottomGrad2 = (fixedPoints.at(5) - fixedPoints.at(3)) / (fixedPoints.at(2) - point2);
                double bottomIntercept2 = fixedPoints.at(5) - (upperGrad2*fixedPoints.at(2));
                double bottomGrad3 = (fixedPoints.at(5) - fixedPoints.at(3)) / (fixedPoints.at(2) - point3);
                double bottomIntercept3 = fixedPoints.at(5) - (upperGrad3*fixedPoints.at(2));
                double rightHandGrad;
                if(fixedPoints.at(1) != fixedPoints.at(2)){
                    rightHandGrad = (fixedPoints.at(4) - fixedPoints.at(5)) / (fixedPoints.at(1) - fixedPoints.at(2));
                }
                else{
                    rightHandGrad = 0;
                }
                double rightHandIntercept = fixedPoints.at(4) - (rightHandGrad*fixedPoints.at(1));
                if(yBinCenter > ((bottomGrad1*xBinCenter) + bottomIntercept1) and yBinCenter < ((upperGrad1*xBinCenter) + upperIntercept1) and xBinCenter >= point1 and ((yBinCenter <= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad < 0) or (yBinCenter >= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad > 0) or (rightHandGrad == 0))){
                    if(signal == "reemitted"){
                        countReEmitted1 = countReEmitted1 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal1 = countTotal1 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter > ((bottomGrad2*xBinCenter) + bottomIntercept2) and yBinCenter < ((upperGrad2*xBinCenter) + upperIntercept2) and xBinCenter >= point2 and ((yBinCenter <= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad < 0) or (yBinCenter >= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad > 0) or (rightHandGrad == 0))){
                    if(signal == "reemitted"){
                        countReEmitted2 = countReEmitted2 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal2 = countTotal2 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter > ((bottomGrad3*xBinCenter) + bottomIntercept3) and yBinCenter < ((upperGrad2*xBinCenter) + upperIntercept2) and xBinCenter >= point3 and ((yBinCenter <= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad < 0) or (yBinCenter >= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad > 0) or (rightHandGrad == 0))){
                    if(signal == "reemitted"){
                        countReEmitted3 = countReEmitted3 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal3 = countTotal3 + allPathsHist->GetBinContent(x,y);                                
                }
            }
            if(numVar == 1){
                double upperGrad1 = (fixedPoints.at(4) - fixedPoints.at(3)) / (point1 - fixedPoints.at(0));
                double upperIntercept1 = fixedPoints.at(4) - (upperGrad1*point1);
                double upperGrad2 = (fixedPoints.at(4) - fixedPoints.at(3)) / (point2 - fixedPoints.at(0));
                double upperIntercept2 = fixedPoints.at(4) - (upperGrad2*point2);
                double upperGrad3 = (fixedPoints.at(4) - fixedPoints.at(3)) / (point3 - fixedPoints.at(0));
                double upperIntercept3 = fixedPoints.at(4) - (upperGrad3*point3);
                double bottomGrad = (fixedPoints.at(5) - fixedPoints.at(3)) / (fixedPoints.at(2) - fixedPoints.at(0));
                double bottomIntercept = fixedPoints.at(5) - (upperGrad1*fixedPoints.at(2));
                double rightHandGrad1;
                double rightHandGrad2;
                double rightHandGrad3;
                if(point1 != fixedPoints.at(2)){
                    rightHandGrad1 = (fixedPoints.at(4) - fixedPoints.at(5)) / (point1 - fixedPoints.at(2));
                }
                else{
                    rightHandGrad1 = 0;
                }
                if(point2 != fixedPoints.at(2)){
                    rightHandGrad2 = (fixedPoints.at(4) - fixedPoints.at(5)) / (point2 - fixedPoints.at(2));
                }
                else{
                    rightHandGrad2 = 0;
                }
                if(point3 != fixedPoints.at(2)){
                    rightHandGrad3 = (fixedPoints.at(4) - fixedPoints.at(5)) / (point3 - fixedPoints.at(2));
                }
                else{
                    rightHandGrad3 = 0;
                }
                double rightHandIntercept1 = fixedPoints.at(5) - (rightHandGrad1*fixedPoints.at(2));
                double rightHandIntercept2 = fixedPoints.at(5) - (rightHandGrad2*fixedPoints.at(2));
                double rightHandIntercept3 = fixedPoints.at(5) - (rightHandGrad3*fixedPoints.at(2));
                if(yBinCenter > ((bottomGrad*xBinCenter) + bottomIntercept) and yBinCenter < ((upperGrad1*xBinCenter) + upperIntercept1) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad1*xBinCenter) + rightHandIntercept1) and rightHandGrad1 < 0) or (yBinCenter >= ((rightHandGrad1*xBinCenter) + rightHandIntercept1) and rightHandGrad1 > 0) or (rightHandGrad1 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted1 = countReEmitted1 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal1 = countTotal1 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter > ((bottomGrad*xBinCenter) + bottomIntercept) and yBinCenter < ((upperGrad2*xBinCenter) + upperIntercept2) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad2*xBinCenter) + rightHandIntercept2) and rightHandGrad2 < 0) or (yBinCenter >= ((rightHandGrad2*xBinCenter) + rightHandIntercept2) and rightHandGrad2 > 0) or (rightHandGrad2 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted2 = countReEmitted2 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal2 = countTotal2 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter > ((bottomGrad*xBinCenter) + bottomIntercept) and yBinCenter < ((upperGrad3*xBinCenter) + upperIntercept3) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad3*xBinCenter) + rightHandIntercept3) and rightHandGrad3 < 0) or (yBinCenter >= ((rightHandGrad3*xBinCenter) + rightHandIntercept3) and rightHandGrad3 > 0) or (rightHandGrad3 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted3 = countReEmitted3 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal3 = countTotal3 + allPathsHist->GetBinContent(x,y);                                
                }
            }
            if(numVar == 2){
                double upperGrad = (fixedPoints.at(4) - fixedPoints.at(3)) / (fixedPoints.at(1) - fixedPoints.at(0));
                double upperIntercept = fixedPoints.at(4) - (upperGrad*fixedPoints.at(1));
                double bottomGrad1 = (fixedPoints.at(5) - fixedPoints.at(3)) / (point1 - fixedPoints.at(0));
                double bottomIntercept1 = fixedPoints.at(5) - (bottomGrad1*point1);
                double bottomGrad2 = (fixedPoints.at(5) - fixedPoints.at(3)) / (point2 - fixedPoints.at(0));
                double bottomIntercept2 = fixedPoints.at(5) - (bottomGrad2*point2);
                double bottomGrad3 = (fixedPoints.at(5) - fixedPoints.at(3)) / (point3 - fixedPoints.at(0));
                double bottomIntercept3 = fixedPoints.at(5) - (bottomGrad3*point3);
                double rightHandGrad1;
                double rightHandGrad2;
                double rightHandGrad3;
                if(point1 != fixedPoints.at(1)){
                    rightHandGrad1 = (fixedPoints.at(4) - fixedPoints.at(5)) / (fixedPoints.at(1) - point1);
                }
                else{
                    rightHandGrad1 = 0;
                }
                if(point2 != fixedPoints.at(1)){
                    rightHandGrad2 = (fixedPoints.at(4) - fixedPoints.at(5)) / (fixedPoints.at(1) - point2);
                }
                else{
                    rightHandGrad2 = 0;
                }
                if(point3 != fixedPoints.at(1)){
                    rightHandGrad3 = (fixedPoints.at(4) - fixedPoints.at(5)) / (fixedPoints.at(1) - point3);
                }
                else{
                    rightHandGrad3 = 0;
                }
                double rightHandIntercept1 = fixedPoints.at(4) - (rightHandGrad1*fixedPoints.at(1));
                double rightHandIntercept2 = fixedPoints.at(4) - (rightHandGrad2*fixedPoints.at(1));
                double rightHandIntercept3 = fixedPoints.at(4) - (rightHandGrad3*fixedPoints.at(1));
                if(yBinCenter > ((bottomGrad1*xBinCenter) + bottomIntercept1) and yBinCenter < ((upperGrad*xBinCenter) + upperIntercept) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad1*xBinCenter) + rightHandIntercept1) and rightHandGrad1 < 0) or (yBinCenter >= ((rightHandGrad1*xBinCenter) + rightHandIntercept1) and rightHandGrad1 > 0) or (rightHandGrad1 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted1 = countReEmitted1 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal1 = countTotal1 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter > ((bottomGrad2*xBinCenter) + bottomIntercept2) and yBinCenter < ((upperGrad*xBinCenter) + upperIntercept) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad2*xBinCenter) + rightHandIntercept2) and rightHandGrad2 < 0) or (yBinCenter >= ((rightHandGrad2*xBinCenter) + rightHandIntercept2) and rightHandGrad2 > 0) or (rightHandGrad2 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted2 = countReEmitted2 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal2 = countTotal2 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter > ((bottomGrad3*xBinCenter) + bottomIntercept3) and yBinCenter < ((upperGrad*xBinCenter) + upperIntercept) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad3*xBinCenter) + rightHandIntercept3) and rightHandGrad3 < 0) or (yBinCenter >= ((rightHandGrad3*xBinCenter) + rightHandIntercept3) and rightHandGrad3 > 0) or (rightHandGrad3 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted3 = countReEmitted3 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal3 = countTotal3 + allPathsHist->GetBinContent(x,y);                                
                }
            }
            if(numVar == 3){
                double bottomGrad1 = (fixedPoints.at(5) - point1) / (fixedPoints.at(2) - fixedPoints.at(0));
                double bottomIntercept1 = fixedPoints.at(5) - (bottomGrad1*fixedPoints.at(2));
                double bottomGrad2 = (fixedPoints.at(5) - point2) / (fixedPoints.at(2) - fixedPoints.at(0));
                double bottomIntercept2 = fixedPoints.at(5) - (bottomGrad2*fixedPoints.at(2));
                double bottomGrad3 = (fixedPoints.at(5) - point3) / (fixedPoints.at(2) - fixedPoints.at(0));
                double bottomIntercept3 = fixedPoints.at(5) - (bottomGrad3*fixedPoints.at(2));
                double upperGrad1 = (fixedPoints.at(4) - point1) / (fixedPoints.at(1) - fixedPoints.at(0));
                double upperIntercept1 = fixedPoints.at(4) - (upperGrad1*fixedPoints.at(1));
                double upperGrad2 = (fixedPoints.at(4) - point2) / (fixedPoints.at(1) - fixedPoints.at(0));
                double upperIntercept2 = fixedPoints.at(4) - (upperGrad2*fixedPoints.at(1));
                double upperGrad3 = (fixedPoints.at(4) - point3) / (fixedPoints.at(1) - fixedPoints.at(0));
                double upperIntercept3 = fixedPoints.at(4) - (upperGrad3*fixedPoints.at(1));
                double rightHandGrad;
                if(fixedPoints.at(1) != fixedPoints.at(2)){
                    rightHandGrad = (fixedPoints.at(4) - fixedPoints.at(5)) / (fixedPoints.at(1) - fixedPoints.at(2));
                }
                else{
                    rightHandGrad = 0;
                }
                double rightHandIntercept = fixedPoints.at(4) - (rightHandGrad*fixedPoints.at(1));
                if(yBinCenter > ((bottomGrad1*xBinCenter) + bottomIntercept1) and yBinCenter < ((upperGrad1*xBinCenter) + upperIntercept1) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad < 0) or (yBinCenter >= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad > 0) or (rightHandGrad == 0))){
                    if(signal == "reemitted"){
                        countReEmitted1 = countReEmitted1 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal1 = countTotal1 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter > ((bottomGrad2*xBinCenter) + bottomIntercept2) and yBinCenter < ((upperGrad2*xBinCenter) + upperIntercept2) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad < 0) or (yBinCenter >= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad > 0) or (rightHandGrad == 0))){
                    if(signal == "reemitted"){
                        countReEmitted2 = countReEmitted2 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal2 = countTotal2 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter > ((bottomGrad3*xBinCenter) + bottomIntercept3) and yBinCenter < ((upperGrad3*xBinCenter) + upperIntercept3) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad < 0) or (yBinCenter >= ((rightHandGrad*xBinCenter) + rightHandIntercept) and rightHandGrad > 0) or (rightHandGrad == 0))){
                    if(signal == "reemitted"){
                        countReEmitted3 = countReEmitted3 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal3 = countTotal3 + allPathsHist->GetBinContent(x,y);                                
                }
            }
            if(numVar == 4){
                double bottomGrad = (fixedPoints.at(5) - fixedPoints.at(3)) / (fixedPoints.at(2) - fixedPoints.at(0));
                double bottomIntercept = fixedPoints.at(5) - (bottomGrad*fixedPoints.at(2));
                double upperGrad1 = (point1 - fixedPoints.at(3)) / (fixedPoints.at(1) - fixedPoints.at(0));
                double upperIntercept1 = point1 - (upperGrad1*fixedPoints.at(1));
                double upperGrad2 = (point2 - fixedPoints.at(3)) / (fixedPoints.at(1) - fixedPoints.at(0));
                double upperIntercept2 = point2 - (upperGrad2*fixedPoints.at(1));
                double upperGrad3 = (point3 - fixedPoints.at(3)) / (fixedPoints.at(1) - fixedPoints.at(0));
                double upperIntercept3 = point3 - (upperGrad3*fixedPoints.at(1));
                double rightHandGrad1;
                double rightHandGrad2;
                double rightHandGrad3;
                if(fixedPoints.at(2) != fixedPoints.at(1)){
                    rightHandGrad1 = (point1 - fixedPoints.at(5)) / (fixedPoints.at(1) - fixedPoints.at(2));
                }
                else{
                    rightHandGrad1 = 0;
                }
                if(fixedPoints.at(2) != fixedPoints.at(1)){
                    rightHandGrad2 = (point2 - fixedPoints.at(5)) / (fixedPoints.at(1) - fixedPoints.at(2));
                }
                else{
                    rightHandGrad2 = 0;
                }
                if(fixedPoints.at(2) != fixedPoints.at(1)){
                    rightHandGrad3 = (point3 - fixedPoints.at(5)) / (fixedPoints.at(1) - fixedPoints.at(2));
                }
                else{
                    rightHandGrad3 = 0;
                }
                double rightHandIntercept1 = fixedPoints.at(5) - (rightHandGrad1*fixedPoints.at(2));
                double rightHandIntercept2 = fixedPoints.at(5) - (rightHandGrad2*fixedPoints.at(2));
                double rightHandIntercept3 = fixedPoints.at(5) - (rightHandGrad3*fixedPoints.at(2));
                if(yBinCenter >= ((bottomGrad*xBinCenter) + bottomIntercept) and yBinCenter <= ((upperGrad1*xBinCenter) + upperIntercept1) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad1*xBinCenter) + rightHandIntercept1) and rightHandGrad1 < 0) or (yBinCenter >= ((rightHandGrad1*xBinCenter) + rightHandIntercept1) and rightHandGrad1 > 0) or (rightHandGrad1 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted1 = countReEmitted1 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal1 = countTotal1 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter >= ((bottomGrad*xBinCenter) + bottomIntercept) and yBinCenter <= ((upperGrad2*xBinCenter) + upperIntercept2) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad2*xBinCenter) + rightHandIntercept2) and rightHandGrad2 < 0) or (yBinCenter >= ((rightHandGrad2*xBinCenter) + rightHandIntercept2) and rightHandGrad2 > 0) or (rightHandGrad2 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted2 = countReEmitted2 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal2 = countTotal2 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter >= ((bottomGrad*xBinCenter) + bottomIntercept) and yBinCenter <= ((upperGrad3*xBinCenter) + upperIntercept3) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad3*xBinCenter) + rightHandIntercept3) and rightHandGrad3 < 0) or (yBinCenter >= ((rightHandGrad3*xBinCenter) + rightHandIntercept3) and rightHandGrad3 > 0) or (rightHandGrad3 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted3 = countReEmitted3 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal3 = countTotal3 + allPathsHist->GetBinContent(x,y);                                
                }
            }
            if(numVar == 5){
                double upperGrad = (fixedPoints.at(4) - fixedPoints.at(3)) / (fixedPoints.at(1) - fixedPoints.at(0));
                double upperIntercept = fixedPoints.at(4) - (upperGrad*fixedPoints.at(1));
                double bottomGrad1 = (point1 - fixedPoints.at(3)) / (fixedPoints.at(2) - fixedPoints.at(0));
                double bottomIntercept1 = point1 - (bottomGrad1*fixedPoints.at(2));
                double bottomGrad2 = (point2 - fixedPoints.at(3)) / (fixedPoints.at(2) - fixedPoints.at(0));
                double bottomIntercept2 = point2 - (bottomGrad2*fixedPoints.at(2));
                double bottomGrad3 = (point3 - fixedPoints.at(3)) / (fixedPoints.at(2) - fixedPoints.at(0));
                double bottomIntercept3 = point3 - (bottomGrad3*fixedPoints.at(2));
                double rightHandGrad1;
                double rightHandGrad2;
                double rightHandGrad3;
                if(fixedPoints.at(2) != fixedPoints.at(1)){
                    rightHandGrad1 = (fixedPoints.at(4) - point1) / (fixedPoints.at(1) - fixedPoints.at(2));
                }
                else{
                    rightHandGrad1 = 0;
                }
                if(fixedPoints.at(2) != fixedPoints.at(1)){
                    rightHandGrad2 = (fixedPoints.at(4) - point2) / (fixedPoints.at(1) - fixedPoints.at(2));
                }
                else{
                    rightHandGrad2 = 0;
                }
                if(fixedPoints.at(2) != fixedPoints.at(1)){
                    rightHandGrad3 = (fixedPoints.at(4) - point3) / (fixedPoints.at(1) - fixedPoints.at(2));
                }
                else{
                    rightHandGrad3 = 0;
                }
                double rightHandIntercept1 = fixedPoints.at(4) - (rightHandGrad1*fixedPoints.at(1));
                double rightHandIntercept2 = fixedPoints.at(4) - (rightHandGrad2*fixedPoints.at(1));
                double rightHandIntercept3 = fixedPoints.at(4) - (rightHandGrad3*fixedPoints.at(1));
                if(yBinCenter > ((bottomGrad1*xBinCenter) + bottomIntercept1) and yBinCenter < ((upperGrad*xBinCenter) + upperIntercept) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad1*xBinCenter) + rightHandIntercept1) and rightHandGrad1 < 0) or (yBinCenter >= ((rightHandGrad1*xBinCenter) + rightHandIntercept1) and rightHandGrad1 > 0) or (rightHandGrad1 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted1 = countReEmitted1 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted1 = countReEmitted1 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal1 = countTotal1 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter > ((bottomGrad2*xBinCenter) + bottomIntercept2) and yBinCenter < ((upperGrad*xBinCenter) + upperIntercept) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad2*xBinCenter) + rightHandIntercept2) and rightHandGrad2 < 0) or (yBinCenter >= ((rightHandGrad2*xBinCenter) + rightHandIntercept2) and rightHandGrad2 > 0) or (rightHandGrad2 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted2 = countReEmitted2 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted2 = countReEmitted2 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal2 = countTotal2 + allPathsHist->GetBinContent(x,y);
                }
                if(yBinCenter > ((bottomGrad3*xBinCenter) + bottomIntercept3) and yBinCenter < ((upperGrad*xBinCenter) + upperIntercept) and xBinCenter >= fixedPoints.at(0) and ((yBinCenter <= ((rightHandGrad3*xBinCenter) + rightHandIntercept3) and rightHandGrad3 < 0) or (yBinCenter >= ((rightHandGrad3*xBinCenter) + rightHandIntercept3) and rightHandGrad3 > 0) or (rightHandGrad3 == 0))){
                    if(signal == "reemitted"){
                        countReEmitted3 = countReEmitted3 + reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted3 = countReEmitted3 + scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal3 = countTotal3 + allPathsHist->GetBinContent(x,y);                                
                }
            }
        }
    }

    double signal_ratio1 = 0;
    double signal_ratio2 = 0;
    double signal_ratio3 = 0;
    if(countReEmitted1 > 0 and countTotal1 > 0) signal_ratio1 = countReEmitted1 / std::sqrt(countTotal1);
    if(countReEmitted2 > 0 and countTotal2 > 0) signal_ratio2 = countReEmitted2 / std::sqrt(countTotal2);
    if(countReEmitted3 > 0 and countTotal3 > 0) signal_ratio3 = countReEmitted3 / std::sqrt(countTotal3);

    std::vector<double> outputFOMs = {signal_ratio1, signal_ratio2, signal_ratio3};

    return outputFOMs;
}

std::vector<double> GetBestFOM(std::vector<double> FOMs, std::vector<double> points){
    std::vector<double> output;
    if(FOMs.at(0) > FOMs.at(1) and FOMs.at(0) > FOMs.at(2)){ //left side best point
        output.push_back(points.at(0));
        if(FOMs.at(1) > FOMs.at(2)){ //right point worst
            output.push_back(points.at(2));
            output.push_back(FOMs.at(0));
            output.push_back(FOMs.at(2));
        }
        else{ // middle point worst
            output.push_back(points.at(1));
            output.push_back(FOMs.at(0));
            output.push_back(FOMs.at(1));
        }
    }
    else if(FOMs.at(1) > FOMs.at(0) and FOMs.at(1) > FOMs.at(2)){ //middle best point
        output.push_back(points.at(1));
        if(FOMs.at(0) > FOMs.at(2)){ //right point worst
            output.push_back(points.at(2));
            output.push_back(FOMs.at(1));
            output.push_back(FOMs.at(2));
        }
        else{ // left point worst
            output.push_back(points.at(0));
            output.push_back(FOMs.at(1));
            output.push_back(FOMs.at(0));
        }
    }
    else if(FOMs.at(2) > FOMs.at(1) and FOMs.at(2) > FOMs.at(0)){ //right side best point
        output.push_back(points.at(2));
        if(FOMs.at(1) > FOMs.at(0)){ //left point worst
            output.push_back(points.at(0));
            output.push_back(FOMs.at(2));
            output.push_back(FOMs.at(0));
        }
        else{ // middle point worst
            output.push_back(points.at(1));
            output.push_back(FOMs.at(2));
            output.push_back(FOMs.at(1));
        }
    }
    else{
        std::cout << "ERROR! No best fit point, points: ";
        for(int i = 0; i< points.size();i++){
            std::cout << points.at(i) << ", ";
        }
        std::cout << std::endl << "FOMs: ";
        for(int i = 0; i< FOMs.size();i++){
            std::cout << FOMs.at(i) << ", ";
        }
    }
    return output;
}

std::vector<TH2F*> GetRegionSelectedHists(std::vector<double> finalPoints, TH2F *hReEmittedPaths, TH2F *hAllPaths, TH2F *hNoisePaths, TH2F *hSingleScatterPaths, TH2F *hOtherPaths, TH2F *hNoEffectPaths, TH2F *hNearReflectPaths, TH2F *hRopesPaths, TH2F *hPMTReflectionPaths, TH2F *hExtWaterScatterPaths, TH2F *hInnerAvReflectPaths, TH2F *hMultipleEffectPaths, TH2F *hAVPipesPaths, TH2F *hAcrylicPaths, TH2F *hOtherScatterPaths){

    //FIXME: pass in vector of hists?

    TH2F *hRegionSelectedReEmittedPaths = (TH2F*)hReEmittedPaths->Clone();
    TH2F *hRegionSelectedAllPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hRegionSelectedNoisePaths = (TH2F*)hNoisePaths->Clone();
    TH2F *hRegionSelectedSingleScatterPaths = (TH2F*)hSingleScatterPaths->Clone();
    TH2F *hRegionSelectedOtherPaths = (TH2F*)hOtherPaths->Clone();
    TH2F *hRegionSelectedNoEffectPaths = (TH2F*)hNoEffectPaths->Clone();
    TH2F *hRegionSelectedNearReflectPaths = (TH2F*)hNearReflectPaths->Clone();
    TH2F *hRegionSelectedRopesPaths = (TH2F*)hRopesPaths->Clone();
    TH2F *hRegionSelectedPMTReflectionPaths = (TH2F*)hPMTReflectionPaths->Clone();
    TH2F *hRegionSelectedExtWaterScatterPaths = (TH2F*)hExtWaterScatterPaths->Clone();
    TH2F *hRegionSelectedInnerAvReflectPaths = (TH2F*)hInnerAvReflectPaths->Clone();
    TH2F *hRegionSelectedMultipleEffectPaths = (TH2F*)hMultipleEffectPaths->Clone();
    TH2F *hRegionSelectedAVPipesPaths = (TH2F*)hAVPipesPaths->Clone();
    TH2F *hRegionSelectedAcrylicPaths = (TH2F*)hAcrylicPaths->Clone();
    TH2F *hRegionSelectedOtherScatterPaths = (TH2F*)hOtherScatterPaths->Clone();

    hRegionSelectedReEmittedPaths->SetName("hRegionSelectedReEmittedPaths");
    hRegionSelectedAllPaths->SetName("hRegionSelectedAllPaths");
    hRegionSelectedNoisePaths->SetName("hRegionCutNoisePaths");
    hRegionSelectedSingleScatterPaths->SetName("hRegionCutSingleScatterPaths");
    hRegionSelectedOtherPaths->SetName("hRegionCutOtherPaths");
    hRegionSelectedNoEffectPaths->SetName("hRegionCutNoEffectPaths");
    hRegionSelectedNearReflectPaths->SetName("hRegionCutNearReflectPaths");
    hRegionSelectedRopesPaths->SetName("hRegionCutRopesPaths");
    hRegionSelectedPMTReflectionPaths->SetName("hRegionCutPMTReflectionPaths");
    hRegionSelectedExtWaterScatterPaths->SetName("hRegionCutExtWaterScatterPaths");
    hRegionSelectedInnerAvReflectPaths->SetName("hRegionCutInnerAvReflectPaths");
    hRegionSelectedMultipleEffectPaths->SetName("hRegionCutMultipleEffectPaths");
    hRegionSelectedAVPipesPaths->SetName("hRegionCutAVPipesPaths");
    hRegionSelectedAcrylicPaths->SetName("hRegionCutAcrylicPaths");
    hRegionSelectedOtherScatterPaths->SetName("hRegionCutOtherScatterPaths");

    TH2F *hDirectCutReEmittedPaths = (TH2F*)hReEmittedPaths->Clone();
    TH2F *hDirectCutAllPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hDirectCutNoisePaths = (TH2F*)hNoisePaths->Clone();
    TH2F *hDirectCutSingleScatterPaths = (TH2F*)hSingleScatterPaths->Clone();
    TH2F *hDirectCutOtherPaths = (TH2F*)hOtherPaths->Clone();
    TH2F *hDirectCutNoEffectPaths = (TH2F*)hNoEffectPaths->Clone();
    TH2F *hDirectCutNearReflectPaths = (TH2F*)hNearReflectPaths->Clone();
    TH2F *hDirectCutRopesPaths = (TH2F*)hRopesPaths->Clone();
    TH2F *hDirectCutPMTReflectionPaths = (TH2F*)hPMTReflectionPaths->Clone();
    TH2F *hDirectCutExtWaterScatterPaths = (TH2F*)hExtWaterScatterPaths->Clone();
    TH2F *hDirectCutInnerAvReflectPaths = (TH2F*)hInnerAvReflectPaths->Clone();
    TH2F *hDirectCutMultipleEffectPaths = (TH2F*)hMultipleEffectPaths->Clone();
    TH2F *hDirectCutAVPipesPaths = (TH2F*)hAVPipesPaths->Clone();
    TH2F *hDirectCutAcrylicPaths = (TH2F*)hAcrylicPaths->Clone();
    TH2F *hDirectCutOtherScatterPaths = (TH2F*)hOtherScatterPaths->Clone();

    hDirectCutReEmittedPaths->SetName("hDirectCutReEmittedPaths");
    hDirectCutAllPaths->SetName("hDirectCutAllPaths");
    hDirectCutNoisePaths->SetName("hDirectCutNoisePaths");
    hDirectCutSingleScatterPaths->SetName("hDirectCutSingleScatterPaths");
    hDirectCutOtherPaths->SetName("hDirectCutOtherPaths");
    hDirectCutNoEffectPaths->SetName("hDirectCutNoEffectPaths");
    hDirectCutNearReflectPaths->SetName("hDirectCutNearReflectPaths");
    hDirectCutRopesPaths->SetName("hDirectCutRopesPaths");
    hDirectCutPMTReflectionPaths->SetName("hDirectCutPMTReflectionPaths");
    hDirectCutExtWaterScatterPaths->SetName("hDirectCutExtWaterScatterPaths");
    hDirectCutInnerAvReflectPaths->SetName("hDirectCutInnerAvReflectPaths");
    hDirectCutMultipleEffectPaths->SetName("hDirectCutMultipleEffectPaths");
    hDirectCutAVPipesPaths->SetName("hDirectCutAVPipesPaths");
    hDirectCutAcrylicPaths->SetName("hDirectCutAcrylicPaths");
    hDirectCutOtherScatterPaths->SetName("hDirectCutOtherScatterPaths");

    TH2F *hReflectedCutReEmittedPaths = (TH2F*)hReEmittedPaths->Clone();
    TH2F *hReflectedCutAllPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutNoisePaths = (TH2F*)hNoisePaths->Clone();
    TH2F *hReflectedCutSingleScatterPaths = (TH2F*)hSingleScatterPaths->Clone();
    TH2F *hReflectedCutOtherPaths = (TH2F*)hOtherPaths->Clone();
    TH2F *hReflectedCutNoEffectPaths = (TH2F*)hNoEffectPaths->Clone();
    TH2F *hReflectedCutNearReflectPaths = (TH2F*)hNearReflectPaths->Clone();
    TH2F *hReflectedCutRopesPaths = (TH2F*)hRopesPaths->Clone();
    TH2F *hReflectedCutPMTReflectionPaths = (TH2F*)hPMTReflectionPaths->Clone();
    TH2F *hReflectedCutExtWaterScatterPaths = (TH2F*)hExtWaterScatterPaths->Clone();
    TH2F *hReflectedCutInnerAvReflectPaths = (TH2F*)hInnerAvReflectPaths->Clone();
    TH2F *hReflectedCutMultipleEffectPaths = (TH2F*)hMultipleEffectPaths->Clone();
    TH2F *hReflectedCutAVPipesPaths = (TH2F*)hAVPipesPaths->Clone();
    TH2F *hReflectedCutAcrylicPaths = (TH2F*)hAcrylicPaths->Clone();
    TH2F *hReflectedCutOtherScatterPaths = (TH2F*)hOtherScatterPaths->Clone();

    hReflectedCutReEmittedPaths->SetName("hReflectedCutReEmittedPaths");
    hReflectedCutAllPaths->SetName("hReflectedCutAllPaths");
    hReflectedCutNoisePaths->SetName("hReflectedCutNoisePaths");
    hReflectedCutSingleScatterPaths->SetName("hReflectedCutSingleScatterPaths");
    hReflectedCutOtherPaths->SetName("hReflectedCutOtherPaths");
    hReflectedCutNoEffectPaths->SetName("hReflectedCutNoEffectPaths");
    hReflectedCutNearReflectPaths->SetName("hReflectedCutNearReflectPaths");
    hReflectedCutRopesPaths->SetName("hReflectedCutRopesPaths");
    hReflectedCutPMTReflectionPaths->SetName("hReflectedCutPMTReflectionPaths");
    hReflectedCutExtWaterScatterPaths->SetName("hReflectedCutExtWaterScatterPaths");
    hReflectedCutInnerAvReflectPaths->SetName("hReflectedCutInnerAvReflectPaths");
    hReflectedCutMultipleEffectPaths->SetName("hReflectedCutMultipleEffectPaths");
    hReflectedCutAVPipesPaths->SetName("hReflectedCutAVPipesPaths");
    hReflectedCutAcrylicPaths->SetName("hReflectedCutAcrylicPaths");
    hReflectedCutOtherScatterPaths->SetName("hReflectedCutOtherScatterPaths");

    double direct_max_time = hNoEffectPaths->ProjectionY()->GetXaxis()->GetBinCenter(hNoEffectPaths->ProjectionY()->GetMaximumBin()) + 10;
    double direct_min_time = hNoEffectPaths->ProjectionY()->GetXaxis()->GetBinCenter(hNoEffectPaths->ProjectionY()->GetMaximumBin()) - 10;
    double direct_cos_alpha = -0.9; //FIXME: don't hardcode this

    double reflected_max_time = hNearReflectPaths->ProjectionY()->GetXaxis()->GetBinCenter(hNearReflectPaths->ProjectionY()->GetMaximumBin()) + 10;
    double reflected_min_time = hNearReflectPaths->ProjectionY()->GetXaxis()->GetBinCenter(hNearReflectPaths->ProjectionY()->GetMaximumBin()) - 10;
    double reflected_cos_alpha = 0.95; //FIXME: don't hardcode this

    std::cout << "Min direct time: " << direct_min_time << std::endl;
    std::cout << "Max direct time: " << direct_max_time << std::endl;
    std::cout << "Min reflected time: " << reflected_min_time << std::endl;
    std::cout << "Max reflected time: " << reflected_max_time << std::endl;

    for(int x=1; x<hReEmittedPaths->GetNbinsX()+1; x++){ //loop over histogram bins
        double xBinCenter = hReEmittedPaths->GetXaxis()->GetBinCenter(x);
        for(int y=1; y<hReEmittedPaths->GetNbinsY()+1; y++){ //loop over histogram bins
            double yBinCenter = hReEmittedPaths->GetYaxis()->GetBinCenter(y);
            double upperGrad = (finalPoints.at(4) - finalPoints.at(3)) / (finalPoints.at(1) - finalPoints.at(0));
            double upperIntercept = finalPoints.at(4) - (upperGrad*finalPoints.at(1));
            double bottomGrad = (finalPoints.at(5) - finalPoints.at(3)) / (finalPoints.at(2) - finalPoints.at(0));
            double bottomIntercept = finalPoints.at(5) - (bottomGrad*finalPoints.at(2));
            double rightHandGrad = (finalPoints.at(4) - finalPoints.at(5)) / (finalPoints.at(1) - finalPoints.at(2));
            double rightHandIntercept = finalPoints.at(4) - (rightHandGrad*finalPoints.at(1));
            if(yBinCenter < ((bottomGrad*xBinCenter) + bottomIntercept) or yBinCenter > ((upperGrad*xBinCenter) + upperIntercept) or xBinCenter < finalPoints.at(0) or (rightHandGrad < 0 and yBinCenter > (rightHandGrad*xBinCenter) + rightHandIntercept) or (rightHandGrad > 0 and yBinCenter < (rightHandGrad*xBinCenter) + rightHandIntercept) or (rightHandGrad == 0 and xBinCenter > finalPoints.at(2))){
                hRegionSelectedReEmittedPaths->SetBinContent(x,y,0);
                hRegionSelectedAllPaths->SetBinContent(x,y,0);
                hRegionSelectedNoisePaths->SetBinContent(x,y,0);
                hRegionSelectedSingleScatterPaths->SetBinContent(x,y,0);
                hRegionSelectedOtherPaths->SetBinContent(x,y,0);
                hRegionSelectedNoEffectPaths->SetBinContent(x,y,0);
                hRegionSelectedNearReflectPaths->SetBinContent(x,y,0);
                hRegionSelectedRopesPaths->SetBinContent(x,y,0);
                hRegionSelectedPMTReflectionPaths->SetBinContent(x,y,0);
                hRegionSelectedExtWaterScatterPaths->SetBinContent(x,y,0);
                hRegionSelectedInnerAvReflectPaths->SetBinContent(x,y,0);
                hRegionSelectedMultipleEffectPaths->SetBinContent(x,y,0);
                hRegionSelectedAVPipesPaths->SetBinContent(x,y,0);
                hRegionSelectedAcrylicPaths->SetBinContent(x,y,0);
                hRegionSelectedOtherScatterPaths->SetBinContent(x,y,0);
            }
            if(xBinCenter > direct_cos_alpha or yBinCenter >= direct_max_time or yBinCenter <= direct_min_time){
                hDirectCutReEmittedPaths->SetBinContent(x,y,0);
                hDirectCutAllPaths->SetBinContent(x,y,0);
                hDirectCutNoisePaths->SetBinContent(x,y,0);
                hDirectCutSingleScatterPaths->SetBinContent(x,y,0);
                hDirectCutOtherPaths->SetBinContent(x,y,0);
                hDirectCutNoEffectPaths->SetBinContent(x,y,0);
                hDirectCutNearReflectPaths->SetBinContent(x,y,0);
                hDirectCutRopesPaths->SetBinContent(x,y,0);
                hDirectCutPMTReflectionPaths->SetBinContent(x,y,0);
                hDirectCutExtWaterScatterPaths->SetBinContent(x,y,0);
                hDirectCutInnerAvReflectPaths->SetBinContent(x,y,0);
                hDirectCutMultipleEffectPaths->SetBinContent(x,y,0);
                hDirectCutAVPipesPaths->SetBinContent(x,y,0);
                hDirectCutAcrylicPaths->SetBinContent(x,y,0);
                hDirectCutOtherScatterPaths->SetBinContent(x,y,0);
            }
            if(xBinCenter < reflected_cos_alpha or yBinCenter >= reflected_max_time or yBinCenter <= reflected_min_time){
                hReflectedCutReEmittedPaths->SetBinContent(x,y,0);
                hReflectedCutAllPaths->SetBinContent(x,y,0);
                hReflectedCutNoisePaths->SetBinContent(x,y,0);
                hReflectedCutSingleScatterPaths->SetBinContent(x,y,0);
                hReflectedCutOtherPaths->SetBinContent(x,y,0);
                hReflectedCutNoEffectPaths->SetBinContent(x,y,0);
                hReflectedCutNearReflectPaths->SetBinContent(x,y,0);
                hReflectedCutRopesPaths->SetBinContent(x,y,0);
                hReflectedCutPMTReflectionPaths->SetBinContent(x,y,0);
                hReflectedCutExtWaterScatterPaths->SetBinContent(x,y,0);
                hReflectedCutInnerAvReflectPaths->SetBinContent(x,y,0);
                hReflectedCutMultipleEffectPaths->SetBinContent(x,y,0);
                hReflectedCutAVPipesPaths->SetBinContent(x,y,0);
                hReflectedCutAcrylicPaths->SetBinContent(x,y,0);
                hReflectedCutOtherScatterPaths->SetBinContent(x,y,0);
            }
        }
    }
    std::vector<TH2F*> outputHists;
    outputHists.push_back(hRegionSelectedReEmittedPaths);
    outputHists.push_back(hRegionSelectedAllPaths);
    outputHists.push_back(hRegionSelectedNoisePaths);
    outputHists.push_back(hRegionSelectedSingleScatterPaths);
    outputHists.push_back(hRegionSelectedOtherPaths);
    outputHists.push_back(hRegionSelectedNoEffectPaths);
    outputHists.push_back(hRegionSelectedNearReflectPaths);
    outputHists.push_back(hRegionSelectedRopesPaths);
    outputHists.push_back(hRegionSelectedPMTReflectionPaths);
    outputHists.push_back(hRegionSelectedExtWaterScatterPaths);
    outputHists.push_back(hRegionSelectedInnerAvReflectPaths);
    outputHists.push_back(hRegionSelectedMultipleEffectPaths);
    outputHists.push_back(hRegionSelectedAVPipesPaths);
    outputHists.push_back(hRegionSelectedAcrylicPaths);
    outputHists.push_back(hRegionSelectedOtherScatterPaths);

    outputHists.push_back(hDirectCutReEmittedPaths);
    outputHists.push_back(hDirectCutAllPaths);
    outputHists.push_back(hDirectCutNoisePaths);
    outputHists.push_back(hDirectCutSingleScatterPaths);
    outputHists.push_back(hDirectCutOtherPaths);
    outputHists.push_back(hDirectCutNoEffectPaths);
    outputHists.push_back(hDirectCutNearReflectPaths);
    outputHists.push_back(hDirectCutRopesPaths);
    outputHists.push_back(hDirectCutPMTReflectionPaths);
    outputHists.push_back(hDirectCutExtWaterScatterPaths);
    outputHists.push_back(hDirectCutInnerAvReflectPaths);
    outputHists.push_back(hDirectCutMultipleEffectPaths);
    outputHists.push_back(hDirectCutAVPipesPaths);
    outputHists.push_back(hDirectCutAcrylicPaths);
    outputHists.push_back(hDirectCutOtherScatterPaths);

    outputHists.push_back(hReflectedCutReEmittedPaths);
    outputHists.push_back(hReflectedCutAllPaths);
    outputHists.push_back(hReflectedCutNoisePaths);
    outputHists.push_back(hReflectedCutSingleScatterPaths);
    outputHists.push_back(hReflectedCutOtherPaths);
    outputHists.push_back(hReflectedCutNoEffectPaths);
    outputHists.push_back(hReflectedCutNearReflectPaths);
    outputHists.push_back(hReflectedCutRopesPaths);
    outputHists.push_back(hReflectedCutPMTReflectionPaths);
    outputHists.push_back(hReflectedCutExtWaterScatterPaths);
    outputHists.push_back(hReflectedCutInnerAvReflectPaths);
    outputHists.push_back(hReflectedCutMultipleEffectPaths);
    outputHists.push_back(hReflectedCutAVPipesPaths);
    outputHists.push_back(hReflectedCutAcrylicPaths);
    outputHists.push_back(hReflectedCutOtherScatterPaths);

    return outputHists;
}

std::vector<double> CheckPoints(std::vector<double> points, std::vector<double> fixedPoints, int numVar){
    if(numVar == 0){
        if(points.at(0) > fixedPoints.at(1)) points.at(0) = fixedPoints.at(1);
        if(points.at(1) > fixedPoints.at(1)) points.at(1) = fixedPoints.at(1);
        if(points.at(2) > fixedPoints.at(1)) points.at(2) = fixedPoints.at(1);
    }
    if(numVar == 1){
        if(points.at(0) < fixedPoints.at(0)) points.at(0) = fixedPoints.at(0);
        if(points.at(1) < fixedPoints.at(0)) points.at(1) = fixedPoints.at(0);
        if(points.at(2) < fixedPoints.at(0)) points.at(2) = fixedPoints.at(0);
    }
    if(numVar == 2){
        if(points.at(0) < fixedPoints.at(0)) points.at(0) = fixedPoints.at(0);
        if(points.at(1) < fixedPoints.at(0)) points.at(1) = fixedPoints.at(0);
        if(points.at(2) < fixedPoints.at(0)) points.at(2) = fixedPoints.at(0);
    }
    if(numVar == 4){
        if(points.at(0) < fixedPoints.at(5)) points.at(0) = fixedPoints.at(5);
        if(points.at(1) < fixedPoints.at(5)) points.at(1) = fixedPoints.at(5);
        if(points.at(2) < fixedPoints.at(5)) points.at(2) = fixedPoints.at(5);
    }
    if(numVar == 5){
        if(points.at(0) > fixedPoints.at(4)) points.at(0) = fixedPoints.at(4);
        if(points.at(1) > fixedPoints.at(4)) points.at(1) = fixedPoints.at(4);
        if(points.at(2) > fixedPoints.at(4)) points.at(2) = fixedPoints.at(4);
    }

    return points;
}