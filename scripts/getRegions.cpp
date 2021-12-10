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
int OptimiseDivideAndConquer(std::string inputFile, int nbins, std::string fibre, bool verbose, bool debug, bool extraInfo, std::string signal_param);
std::vector<double> GetThreePoints(double bestPoint, double worstPoint, std::vector<double> originalPoints);
std::vector<double> GetFOMs(std::vector<double> points, std::vector<double> fixedPoints, int numVar, TH2F *allPathsHist, TH2F *reEmittedHist, TH2F *scatteredHist, std::string signal);
std::vector<double> GetBestFOM(std::vector<double> FOMs, std::vector<double> points);
std::vector<TH2F*> GetRegionSelectedHists(std::vector<double> finalPoints, TH2F *hReEmittedPaths, TH2F *hAllPaths, TH2F *hNoisePaths, TH2F *hSingleScatterPaths, TH2F *hOtherPaths, TH2F *hNoEffectPaths, TH2F *hNearReflectPaths, TH2F *hRopesPaths, TH2F *hPMTReflectionPaths, TH2F *hExtWaterScatterPaths, TH2F *hInnerAvReflectPaths, TH2F *hMultipleEffectPaths, TH2F *hAVPipesPaths, TH2F *hAcrylicPaths, TH2F *hOtherScatterPaths, std::string saveroot_txt);
std::vector<double> CheckPoints(std::vector<double> points, std::vector<double> fixedPoints, int numVar);

int main(int argc, char** argv){
    std::string file = argv[1];
    int nbins = std::stoi(argv[2]);
    std::string fibre = argv[3];
    bool verbose = std::stoi(argv[4]);
    bool debug = std::stoi(argv[5]);
    bool extraInfo = std::stoi(argv[6]);
    std::string signal = argv[7];
    auto t1 = std::chrono::high_resolution_clock::now();
    int status = OptimiseDivideAndConquer(file, nbins, fibre, verbose, debug, extraInfo, signal); 
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "Script took " << std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() / 1E6 << " s to execute" << std::endl;

    return 0;
}

/**
 * @brief Finds Triangular region with best (highest) FOM (attenuated / sqrt(total)) via a divide and conquer method
 * Start with p=(min, mid, max) for each point (x and y coordinates), as well as initial fixed points that are the largest
 * right-angle triangle to start with. Find FOM for triangle by replacing each point with each option in p (one at a time),
 * while keeping the other fixed points temprarily fixed. Then replace the associated fixed point with whichever of min or
 * max gave the best FOM. mid then replaces whichever one was used, and a new mid is calculated. Continue until changes is
 * less than tolerance.
 * 
 * @param inputFile Root file with original histograms with all the data.
 * @param nbins Number of bins (resolution) that we wish to use (max is capped at number of bins of original hists).
 * @param fibre Fibre used to determine its angle offset and thus the direct beam spot angle.
 * @param verbose Print extra info.
 * @param debug Print extra info.
 * @param extraInfo Print extra info.
 * @param signal_param signal = reemitted, scattered or attenuated. The signal we are trying to optimise for.
 * @return int 
 */
int OptimiseDivideAndConquer(std::string inputFile, int nbins, std::string fibre, bool verbose, bool debug, bool extraInfo, std::string signal_param){

    auto timeStart = std::chrono::high_resolution_clock::now();

    if(signal_param != "reemitted" and signal_param != "scattered" and signal_param != "attenuated"){ // attenuated = scattered + reemitted
        std::cout << "Invalid signal! Use: reemitted, scattered, attenuated";
        return 1;
    }

    //create histograms and graphs

    TH2F *hReEmittedPaths; TH2F *hAllPaths; TH2F *hNoisePaths;
    TH2F *hSingleScatterPaths; TH2F *hOtherPaths; TH2F *hNoEffectPaths;
    TH2F *hNearReflectPaths; TH2F *hRopesPaths; TH2F *hPMTReflectionPaths;
    TH2F *hExtWaterScatterPaths; TH2F *hInnerAvReflectPaths; TH2F *hMultipleEffectPaths;
    TH2F *hAVPipesPaths; TH2F *hAcrylicPaths; TH2F *hOtherScatterPaths;

    std::vector<TH2F*> Hists;
    Hists.push_back(hReEmittedPaths); Hists.push_back(hAllPaths);
    Hists.push_back(hNoisePaths); Hists.push_back(hSingleScatterPaths);
    Hists.push_back(hOtherPaths); Hists.push_back(hNoEffectPaths);
    Hists.push_back(hNearReflectPaths); Hists.push_back(hRopesPaths);
    Hists.push_back(hPMTReflectionPaths); Hists.push_back(hExtWaterScatterPaths);
    Hists.push_back(hInnerAvReflectPaths); Hists.push_back(hMultipleEffectPaths);
    Hists.push_back(hAVPipesPaths); Hists.push_back(hAcrylicPaths);
    Hists.push_back(hOtherScatterPaths);

    TGraph *hFOMx_a = new TGraph(); TGraph *hFOMx_b = new TGraph();
    TGraph *hFOMx_c = new TGraph(); TGraph *hFOMy_a = new TGraph();
    TGraph *hFOMy_b = new TGraph(); TGraph *hFOMy_c = new TGraph();
    TGraph *hFOMmainloop = new TGraph(); TGraph *hDiffx_a = new TGraph();
    TGraph *hDiffx_b = new TGraph(); TGraph *hDiffx_c = new TGraph();
    TGraph *hDiffy_a = new TGraph(); TGraph *hDiffy_b = new TGraph();
    TGraph *hDiffy_c = new TGraph(); TGraph *hPointx_a = new TGraph();
    TGraph *hPointx_b = new TGraph(); TGraph *hPointx_c = new TGraph();
    TGraph *hPointy_a = new TGraph(); TGraph *hPointy_b = new TGraph();
    TGraph *hPointy_c = new TGraph();

    std::vector<TGraph*> hFOMxy;
    hFOMxy.push_back(hFOMx_a); hFOMxy.push_back(hFOMx_b);
    hFOMxy.push_back(hFOMx_c); hFOMxy.push_back(hFOMy_a);
    hFOMxy.push_back(hFOMy_b); hFOMxy.push_back(hFOMy_c);

    std::vector<TGraph*> hDiffxy;
    hDiffxy.push_back(hDiffx_a); hDiffxy.push_back(hDiffx_b);
    hDiffxy.push_back(hDiffx_c); hDiffxy.push_back(hDiffy_a);
    hDiffxy.push_back(hDiffy_b); hDiffxy.push_back(hDiffy_c);
    
    std::vector<TGraph*> hPointxy;
    Graphs.push_back(hPointx_a); Graphs.push_back(hPointx_b);
    Graphs.push_back(hPointx_c); Graphs.push_back(hPointy_a);
    Graphs.push_back(hPointy_b); Graphs.push_back(hPointy_c);

    std::vector<std::string> Hist_names;
    Hist_names.push_back("hReemissionResTimeVsCosTheta");
    Hist_names.push_back("hPmtResTimeVsCosTheta");
    Hist_names.push_back("hNoiseResTimeVsCosTheta");
    Hist_names.push_back("hSingleScatterResTimeVsCosTheta");
    Hist_names.push_back("hOtherEffectResTimeVsCosTheta");
    Hist_names.push_back("hNoEffectResTimeVsCosTheta");
    Hist_names.push_back("hNearReflectResTimeVsCosTheta");
    Hist_names.push_back("hRopesResTimeVsCosTheta");
    Hist_names.push_back("hPMTReflectionResTimeVsCosTheta");
    Hist_names.push_back("hExtWaterScatterResTimeVsCosTheta");
    Hist_names.push_back("hInnerAvReflectionResTimeVsCosTheta");
    Hist_names.push_back("hMultipleEffectResTimeVsCosTheta");
    Hist_names.push_back("hAVPipesResTimeVsCosTheta");
    Hist_names.push_back("hAcrylicScatterResTimeVsCosTheta");
    Hist_names.push_back("OtherScatterResTimeVsCosTheta");

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
    for (i = 0; i < Hist_names.size(); ++i) {
        raw_file->GetObject(Hist_names.at(i), Hists.at(i));
    }

    //set up constants
    double countCalculations = 0;
    double countCycles = 0;

    std::vector<std::string> point_names = {"x_a", "x_b", "x_c", "y_a", "y_b", "y_c"};
    double y_min = hReEmittedPaths->GetYaxis()->GetXmin();
    double y_max = hReEmittedPaths->GetYaxis()->GetXmax();
    double point_mins[6] = {-1, -1, -1, y_min, y_min, y_min};
    double point_maxs[6] = {1, 1, 1, y_max, y_max, y_max};
    std::vector<double> point_tolerances;
    for (int i = 0; i < 6; ++i) {
        point_tolerances.push_back((point_maxs[i] - point_mins[i]) / nbins);
    }

    int nBinsX = hReEmittedPaths->GetXaxis()->GetNbins();
    int nBinsY = hReEmittedPaths->GetYaxis()->GetNbins();

    double xBinWidth = hReEmittedPaths->GetXaxis()->GetBinCenter(2) - hReEmittedPaths->GetXaxis()->GetBinCenter(1);
    double yBinWidth = hReEmittedPaths->GetYaxis()->GetBinCenter(2) - hReEmittedPaths->GetYaxis()->GetBinCenter(1);

    for (int i = 0; i < 3; ++i) {
        if (point_tolerances.at(i) < xBinWidth) {
            point_tolerances.at(i) = xBinWidth;
            std::cout << point_names.at(i) << "_tolerance is set to bin width, nbins: " << (point_maxs[i] - point_mins[i]) / xBinWidth << std::endl;
        }
        j = i + 3;
        if (point_tolerances.at(j) < yBinWidth) {
            point_tolerances.at(j) = yBinWidth;
            std::cout << point_names.at(j) << "_tolerance is set to bin width, nbins: " << (point_maxs[j] - point_mins[j]) / yBinWidth << std::endl;
        }
    }

    if(debug) std::cout << "x_a_tolerance: " << point_tolerances.at(0) << ", x_b_tolerance: " << point_tolerances.at(1) << std::endl;

    std::vector<double> points_diffs = {9999999999, 9999999999, 9999999999, 9999999999, 9999999999, 9999999999};
    std::vector<double> points_temp_diffs = {9999999999, 9999999999, 9999999999, 9999999999, 9999999999, 9999999999};

    auto timeSetup = std::chrono::high_resolution_clock::now();

    if(verbose) std::cout << "Time to read in file and set up hists: " << std::chrono::duration_cast<std::chrono::microseconds>(timeSetup - timeStart).count() / 1E6 << "s" << std::endl;
    // do loop: x_a, x_b, x_c, y_a, y_b, y_c
    // points: a vector containing the points to change for three regions
    // fixedPoints: a vector containing all the static points

    auto timeStartLoop = std::chrono::high_resolution_clock::now();

    // Initialise variables:
    // fixedPoints = {x_a_min, x_b_max, x_c_max, y_a_min, y_b_max, y_c_min};
    std::vector<double> fixedPoints = {point_mins[0], point_maxs[1], point_maxs[2], point_mins[3], point_maxs[4], point_mins[5]};
    bool firstRun = true;

    int numMainLoopIterations = 0;
    int num_iterations[6] = {0, 0, 0, 0, 0, 0};
    double prev_best_points_main[6];
    bool first_runs[6];
    double prevBestFOMPoints[6];
    double temp_diffs[6];
    std::vector<double> FOMs[6] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    std::vector<double> bestworstFOMPoints[6] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

    std::vector<double> points[6];
    for (int i = 0; i < 6; ++i) {
        points[i] = {point_mins[i], (point_maxs[i] + point_mins[i])/2., point_maxs[i]};
    }
    
    //while(((x_a_diff > x_a_tolerance or x_b_diff > x_b_tolerance or x_c_diff > x_c_tolerance or y_a_diff > y_a_tolerance or y_b_diff > y_b_tolerance or y_c_diff > y_c_tolerance) and !firstRun) or firstRun){
    while((points_diffs > point_tolerances and !firstRun) or firstRun){
        first_runs[6] = {true, true, true, true, true, true};
        prevBestFOMPoints[6] = {9999999999, 9999999999, 9999999999, 9999999999, 9999999999, 9999999999};
        temp_diffs[6] = {9999999999, 9999999999, 9999999999, 9999999999, 9999999999, 9999999999};

        if(debug) std::cout << "In main while loop with x_a_temp_diff: " << std::abs(temp_diffs[0]) << std::endl;

        for (int i = 0; i < 6; ++i) {
            while(std::abs(temp_diffs[i]) > point_tolerances.at(i)){
                if(first_runs[i]){
                    points[i] = {point_mins[i], (point_maxs[i] + point_mins[i])/2., point_maxs[i]}; //l, c, r
                    points[i] = CheckPoints(points[i], fixedPoints, 0);
                    first_runs[i] = false;
                    if(debug) std::cout << "Set up first " << point_names[i] << "_run points: " << points[i].at(0)
                                        << ", " << points[i].at(1) << ", " << points[i].at(2) << std::endl;
                }
                FOMs.at(i) = GetFOMs(points[i], fixedPoints, i, hAllPaths, hReEmittedPaths, hSingleScatterPaths, signal_param);
                if(debug) std::cout << "Got FOMs " << FOMs.at(i).at(0) << ", " << FOMs.at(i).at(1) << ", " << FOMs.at(i).at(2) << std::endl;
                
                bestworstFOMPoints[i] = GetBestFOM(FOMs.at(i), points[i]);
                if(debug) std::cout << "Got best FOMs" << std::endl;
                
                points[i] = GetThreePoints(bestworstFOMPoints[i].at(0), bestworstFOMPoints[i].at(1), points[i]);
                points[i] = CheckPoints(points[i], fixedPoints, i);
                if(debug) std::cout << "Got new points: " << points[i].at(0) << ", " << points[i].at(1) << ", " << points[i].at(2) << std::endl;
                
                temp_diffs[i] = std::abs(points[i].at(1) - prevBestFOMPoints[i]);
                if(debug) std::cout << "Got difference: " << temp_diffs[i] << std::endl;
                
                prevBestFOMPoints[i] = points[i].at(1);
                if(extraInfo) hFOMxy.at(i)->SetPoint(num_iterations[i], num_iterations[i], bestworstFOMPoints[i].at(2));
                num_iterations[i]++;
            }
        }

        if(debug) std::cout << "Setting fixed point: " << bestworstFOMPoints[5].at(0) << std::endl;
        fixedPoints.at(5) = bestworstFOMPoints[5].at(0);
        if(debug) std::cout << "" << std::endl;

        if(debug) std::cout << "Getting main differences" << std::endl;

        for (int i = 0; i < 6; ++i) {
            if(points_diffs.at(i) == 9999999999){
                points_diffs.at(i) = std::abs(bestworstFOMPoints[i].at(0));
                prev_best_points_main[i] = bestworstFOMPoints[i].at(0);
            }
            else{
                points_diffs.at(i) = std::abs(prev_best_points_main[i] - bestworstFOMPoints[i].at(0));
                prev_best_points_main[i] = bestworstFOMPoints[i].at(0);
                firstRun = false;
            }
        }
        
        if(debug) std::cout << "Got main differences" << std::endl;
        if(debug) std::cout << "" << std::endl;
        if(extraInfo) {
            hFOMmainloop->SetPoint(numMainLoopIterations, numMainLoopIterations, bestworstFOMPoints[5].at(2));
            for (int i = 0; i < 6; ++i) {
                hDiffxy.at(i)->SetPoint(numMainLoopIterations, numMainLoopIterations, points_diffs.at(i));
                hPointxy.at(i)->SetPoint(numMainLoopIterations, numMainLoopIterations, fixedPoints.at(i));
            }
        }
        numMainLoopIterations++;
    }

    // open txt file and print Region limits (can then use on real data!)

    std::ofstream outputFile_txt;
    // get file name from path+filename string
    std::size_t botDirPos_txt = inputFile.find_last_of("/");
    std::string filename_txt = inputFile.substr(botDirPos_txt+1, inputFile.length()) + ".txt";
    std::string saveroot_txt = "region_selected_limits_" + signal_param + "_" + filename_txt;
    outputFile_txt.open(saveroot_txt.c_str());
    outputFile_txt << std::to_string(fixedPoints.at(0)) + "\n"; // x_a
    outputFile_txt << std::to_string(fixedPoints.at(1)) + "\n"; // x_b
    outputFile_txt << std::to_string(fixedPoints.at(2)) + "\n"; // x_c
    outputFile_txt << std::to_string(fixedPoints.at(3)) + "\n"; // y_a
    outputFile_txt << std::to_string(fixedPoints.at(4)) + "\n"; // y_b
    outputFile_txt << std::to_string(fixedPoints.at(5)) + "\n"; // y_c
    outputFile_txt.close();

    // now use final result to get region, write to file

    if(verbose) std::cout << "Region found" << std::endl;
    if(verbose) std::cout << "    x_a: " << fixedPoints.at(0) << std::endl;
    if(verbose) std::cout << "    x_b: " << fixedPoints.at(1) << std::endl;
    if(verbose) std::cout << "    x_c: " << fixedPoints.at(2) << std::endl;
    if(verbose) std::cout << "    y_a: " << fixedPoints.at(3) << std::endl;
    if(verbose) std::cout << "    y_b: " << fixedPoints.at(4) << std::endl;
    if(verbose) std::cout << "    y_c: " << fixedPoints.at(5) << std::endl;
    if(verbose) std::cout << "There were " << numMainLoopIterations << " iterations of the main loop" << std::endl;

    std::vector<TH2F*> regionSelectedHists = GetRegionSelectedHists(fixedPoints, hReEmittedPaths, hAllPaths, hNoisePaths, hSingleScatterPaths, hOtherPaths, hNoEffectPaths, hNearReflectPaths, hRopesPaths, hPMTReflectionPaths, hExtWaterScatterPaths, hInnerAvReflectPaths, hMultipleEffectPaths, hAVPipesPaths, hAcrylicPaths, hOtherScatterPaths, saveroot_txt);

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

    // get points
    double x_a = fixedPoints.at(0);
    double x_b = fixedPoints.at(1);
    double x_c = fixedPoints.at(2);
    double y_a = fixedPoints.at(3);
    double y_b = fixedPoints.at(4);
    double y_c = fixedPoints.at(5);

    double direct_max_time = hNoEffectPaths->ProjectionY()->GetXaxis()->GetBinCenter(hNoEffectPaths->ProjectionY()->GetMaximumBin()) + 10;
    double direct_min_time = hNoEffectPaths->ProjectionY()->GetXaxis()->GetBinCenter(hNoEffectPaths->ProjectionY()->GetMaximumBin()) - 10;
    //FIXME: don't hardcode this:
    double direct_cos_alpha;
    if (fibre == "FA089") {  // 10deg off-axis
        direct_cos_alpha = -0.85; 
    } else if (fibre == "FA173" or fibre == "FA150" or fibre == "FA093") {  // 20deg off-axis
        direct_cos_alpha = -0.6;
    } else {  // on-axis
        direct_cos_alpha = -0.9;
    }

    double reflected_max_time = hNearReflectPaths->ProjectionY()->GetXaxis()->GetBinCenter(hNearReflectPaths->ProjectionY()->GetMaximumBin()) + 10;
    double reflected_min_time = hNearReflectPaths->ProjectionY()->GetXaxis()->GetBinCenter(hNearReflectPaths->ProjectionY()->GetMaximumBin()) - 10;
    double reflected_cos_alpha = 0.95; //FIXME: don't hardcode this

    // Draw box cuts on resthit vs costheta hist
    TCanvas *c1 = new TCanvas("cuts","cuts");  //Create output canvas to be saved in output file
    TH2F *h = (TH2F*)hAllPaths->Clone();
    h->Draw("colz");  // Draw histogram

    // create lines
    std::vector<TLine> lines;
    // region
    lines.push_back(TLine(x_a, y_a, x_b, y_b));
    lines.push_back(TLine(x_a, y_a, x_c, y_c));
    lines.push_back(TLine(x_c, y_c, x_b, y_b));
    // direct box
    lines.push_back(TLine(direct_cos_alpha, direct_min_time, direct_cos_alpha, direct_max_time));
    lines.push_back(TLine(-1, direct_min_time, -1, direct_max_time));
    lines.push_back(TLine(-1, direct_max_time, direct_cos_alpha, direct_max_time));
    lines.push_back(TLine(-1, direct_min_time, direct_cos_alpha, direct_min_time));
    // reflected box
    lines.push_back(TLine(reflected_cos_alpha, reflected_min_time, reflected_cos_alpha, reflected_max_time));
    lines.push_back(TLine(1, reflected_min_time, 1, reflected_max_time));
    lines.push_back(TLine(1, reflected_max_time, reflected_cos_alpha, reflected_max_time));
    lines.push_back(TLine(1, reflected_min_time, reflected_cos_alpha, reflected_min_time));

    // draw lines
    for(int i=0; i<lines.size(); ++i){
        lines[i].SetLineColor(kBlack);
        lines[i].Draw("SAME");
    }

    // Write canvas to root file
    rootfile->cd();
    c1->Write();  
    delete c1;

    rootfile->Write();
    rootfile->Close();

    return 0;
}

/**
 * @brief Remove the point in points with the worst FOM, and recalculate the mid point from the two remaining.
 * If the middle point is the worst, remove the second worst instead.
 * 
 * @param bestPoint Point with best FOM
 * @param worstPoint Point with worst FOM
 * @param originalPoints (z_i_min, z_i_mid, z_i_max), for z=x,y, and i in {a, b, c}
 * @return std::vector<double> 
 */
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


/**
 * @brief Replace the associated point in the trianle (normally given by fixedPoints) with each point (of 3) in points,
 * draw a triangle with each such configuration, and count the number of different types on events that fall inside/outside
 * of the triangle. The FOM for each is a ratio of numbers of events: signal inside / sqrt(total inside), where signal is
 * defined below, and total means any event type.
 * 
 * @param points (z_i_min, z_i_mid, z_i_max), for z=x,y, and i in {a, b, c}
 * @param fixedPoints  (x_a_min, x_b_max, x_c_max, y_a_min, y_b_max, y_c_min)
 * @param numVar Denotes which point is being used/replaced (0=x_a, 1=x_b, 2=x_c, 3=y_a, 4=y_b, 5=y_c)
 * @param allPathsHist Histogram for all light paths
 * @param reEmittedHist Histogram for reemitted light paths
 * @param scatteredHist Histogram for scattered light paths
 * @param signal = scattered, reemitted or attenuated. (attenuated = scattered + reemitted). Note that in the code countReEmitted(1,2,3)
 * refers to whichever one of these was selected.
 * @return std::vector<double> 
 */
std::vector<double> GetFOMs(std::vector<double> points, std::vector<double> fixedPoints, int numVar, TH2F *allPathsHist, TH2F *reEmittedHist, TH2F *scatteredHist, std::string signal){

    //FIXME: Pass in vector of hists?

    double countReEmitted[3] = {0, 0, 0};
    double countTotal[3] = {0, 0, 0};

    // Create triangle with fixed points
    triangle Tri = triangle(fixedPoints.at(0), fixedPoints.at(1), fixedPoints.at(2), fixedPoints.at(3),
                            fixedPoints.at(4), fixedPoints.at(5));

    for(int x=1; x<reEmittedHist->GetNbinsX()+1; x++){ //loop over histogram bins
        double xBinCenter = reEmittedHist->GetXaxis()->GetBinCenter(x);
        for(int y=1; y<reEmittedHist->GetNbinsY()+1; y++){
            double yBinCenter = reEmittedHist->GetYaxis()->GetBinCenter(y);
            // Replace the appropriate point in the triangle with each point in points and see if the bin falls in the triangle.
            for (int i = 0; i < 3; ++i){
                Tri[numVar] = points.at(i);
                if(Tri.check_point_inside_triangle(xBinCenter, yBinCenter)){
                    if(signal == "reemitted"){
                        countReEmitted[i] += reEmittedHist->GetBinContent(x,y);
                    }
                    else if(signal == "scattered"){
                        countReEmitted[i] += scatteredHist->GetBinContent(x,y);
                    }
                    else if(signal == "attenuated"){
                        countReEmitted[i] += scatteredHist->GetBinContent(x,y) + reEmittedHist->GetBinContent(x,y);
                    }
                    countTotal[i] += allPathsHist->GetBinContent(x,y);
                }
            }
        }
    }

    double signal_ratios[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i){
        if(countReEmitted[i] > 0 and countTotal[i] > 0) {
            signal_ratios[i] = countReEmitted[i] / std::sqrt(countTotal[i]);
        }
    }

    std::vector<double> outputFOMs = {signal_ratios[0], signal_ratios[1], signal_ratios[2]};
    return outputFOMs;
}


/**
 * @brief Compares the FOM of each point (of 3) in points, and returns a vector with:
 * (point with best (highest) FOM, point with worst (lowest) FOM, FOM of best point, FOM of worst point).
 * 
 * @param FOMs Figures of merit, See GetFOM() function.
 * @param points (z_i_min, z_i_mid, z_i_max), for z=x,y, and i in {a, b, c}
 * @return std::vector<double> 
 */
std::vector<double> GetBestFOM(std::vector<double> FOMs, std::vector<double> points){
    std::vector<double> output;
    if(FOMs.at(0) >= FOMs.at(1) and FOMs.at(0) >= FOMs.at(2)){ //left side best point
        output.push_back(points.at(0));
        if(FOMs.at(1) >= FOMs.at(2)){ //right point worst
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
    else if(FOMs.at(1) >= FOMs.at(0) and FOMs.at(1) >= FOMs.at(2)){ //middle best point
        output.push_back(points.at(1));
        if(FOMs.at(0) >= FOMs.at(2)){ //right point worst
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
    else if(FOMs.at(2) >= FOMs.at(1) and FOMs.at(2) >= FOMs.at(0)){ //right side best point
        output.push_back(points.at(2));
        if(FOMs.at(1) >= FOMs.at(0)){ //left point worst
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

/**
 * @brief Clone all the original histograms (all but the first params below) 3 times: One for the triangle region,
 * called Region, and one for the direct (called Direct) and reflected (Reflected) beam spots each.
 * The bins falling outside of these regions are then set to zero in each corresponding histogram.
 * Returns vector of resultings histograms.
 * 
 * @param finalPoints final fixed points (x_a_min, x_b_max, x_c_max, y_a_min, y_b_max, y_c_min).
 * @param hReEmittedPaths Original histogram.
 * @param hAllPaths Original histogram.
 * @param hNoisePaths Original histogram.
 * @param hSingleScatterPaths Original histogram.
 * @param hOtherPaths Original histogram.
 * @param hNoEffectPaths Original histogram.
 * @param hNearReflectPaths Original histogram.
 * @param hRopesPaths Original histogram.
 * @param hPMTReflectionPaths Original histogram.
 * @param hExtWaterScatterPaths Original histogram.
 * @param hInnerAvReflectPaths Original histogram.
 * @param hMultipleEffectPaths Original histogram.
 * @param hAVPipesPaths Original histogram.
 * @param hAcrylicPaths Original histogram.
 * @param hOtherScatterPaths Original histogram.
 * @return std::vector<TH2F*> 
 */
std::vector<TH2F*> GetRegionSelectedHists(std::vector<double> finalPoints, TH2F *hReEmittedPaths, TH2F *hAllPaths, TH2F *hNoisePaths, TH2F *hSingleScatterPaths, TH2F *hOtherPaths, TH2F *hNoEffectPaths, TH2F *hNearReflectPaths, TH2F *hRopesPaths, TH2F *hPMTReflectionPaths, TH2F *hExtWaterScatterPaths, TH2F *hInnerAvReflectPaths, TH2F *hMultipleEffectPaths, TH2F *hAVPipesPaths, TH2F *hAcrylicPaths, TH2F *hOtherScatterPaths, std::string saveroot_txt){

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

    // Print to Regions lims file (append)
    std::ofstream outputFile_txt;
    outputFile_txt.open(saveroot_txt.c_str(), std::ios::app);
    outputFile_txt << std::to_string(direct_max_time) + "\n";
    outputFile_txt << std::to_string(direct_min_time) + "\n";
    outputFile_txt << std::to_string(direct_cos_alpha) + "\n";
    outputFile_txt << std::to_string(reflected_max_time) + "\n";
    outputFile_txt << std::to_string(reflected_min_time) + "\n";
    outputFile_txt << std::to_string(reflected_cos_alpha) + "\n";

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


/**
 * @brief Checks if points have broken the triangle (for ex if x_a>x_b, which would flip the order around).
 * If so replace with associated limit from fixedPoints.
 * 
 * @param points (z_i_min, z_i_mid, z_i_max), for z=x,y, and i in {a, b, c}.
 * @param fixedPoints (x_a_min, x_b_max, x_c_max, y_a_min, y_b_max, y_c_min).
 * @param numVar Denotes which point is being checked (0=x_a, 1=x_b, 2=x_c, 3=y_a, 4=y_b, 5=y_c).
 * @return std::vector<double> 
 */
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