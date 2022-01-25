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
void DrawRegionLims(std::vector<double> fixedPoints, HistList hists_lists, std::string fibre);
std::vector<double> GetFOMs(std::vector<double> points, std::vector<double> fixedPoints, int numVar, TH2F *allPathsHist, TH2F *reEmittedHist, TH2F *scatteredHist, std::string signal);
std::vector<double> GetBestFOM(std::vector<double> FOMs, std::vector<double> points);
HistList GetRegionSelectedHists(std::vector<double> finalPoints, HistList hists_lists, std::string fibre, std::string saveroot_txt);

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
    HistList hists_lists = HistList(inputFile);

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
    hPointxy.push_back(hPointx_a); hPointxy.push_back(hPointx_b);
    hPointxy.push_back(hPointx_c); hPointxy.push_back(hPointy_a);
    hPointxy.push_back(hPointy_b); hPointxy.push_back(hPointy_c);

    //set up constants
    std::string point_names[6] = {"x_a", "x_b", "x_c", "y_a", "y_b", "y_c"};
    double y_min = hists_lists.Tracking_Hists().at(0)->GetYaxis()->GetXmin();
    double y_max = hists_lists.Tracking_Hists().at(0)->GetYaxis()->GetXmax();
    double point_mins[6] = {-1, -1, -1, y_min, y_min, y_min};
    double point_maxs[6] = {1, 1, 1, y_max, y_max, y_max};
    std::vector<double> point_tolerances;
    for (int i = 0; i < 6; ++i) {
        point_tolerances.push_back((point_maxs[i] - point_mins[i]) / nbins);
    }

    double xBinWidth = hists_lists.Tracking_Hists().at(0)->GetXaxis()->GetBinCenter(2) - hists_lists.Tracking_Hists().at(0)->GetXaxis()->GetBinCenter(1);
    double yBinWidth = hists_lists.Tracking_Hists().at(0)->GetYaxis()->GetBinCenter(2) - hists_lists.Tracking_Hists().at(0)->GetYaxis()->GetBinCenter(1);

    int j;
    for (int i = 0; i < 3; ++i) {
        if (point_tolerances.at(i) < xBinWidth) {
            point_tolerances.at(i) = xBinWidth;
            std::cout << point_names[i] << "_tolerance is set to bin width, nbins: " << (point_maxs[i] - point_mins[i]) / xBinWidth << std::endl;
        }
        j = i + 3;
        if (point_tolerances.at(j) < yBinWidth) {
            point_tolerances.at(j) = yBinWidth;
            std::cout << point_names[j] << "_tolerance is set to bin width, nbins: " << (point_maxs[j] - point_mins[j]) / yBinWidth << std::endl;
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

    bool loop_condition = true; // If any point difference is still larger than its associated tolerance, keep looping. See end of while loop.
    while(loop_condition or firstRun){
        // Set/Reset variables
        for (int i = 0; i < 6; ++i) {
            first_runs[i] = true;
            prevBestFOMPoints[i] = 9999999999;
            temp_diffs[i] = 9999999999;
        }

        if(debug) std::cout << "In main while loop with x_a_temp_diff: " << std::abs(temp_diffs[0]) << std::endl;

        for (int i = 0; i < 6; ++i) {
            while(std::abs(temp_diffs[i]) > point_tolerances.at(i)){
                if(first_runs[i]){
                    points[i] = {point_mins[i], (point_maxs[i] + point_mins[i])/2., point_maxs[i]}; //l, c, r
                    first_runs[i] = false;
                    if(debug) std::cout << "Set up first " << point_names[i] << "_run points: " << points[i].at(0)
                                        << ", " << points[i].at(1) << ", " << points[i].at(2) << std::endl;
                }
                FOMs[i] = GetFOMs(points[i], fixedPoints, i, hists_lists.Tracking_Hists().at(1), hists_lists.Tracking_Hists().at(0), hists_lists.Tracking_Hists().at(3), signal_param);  //hAllPaths, hReEmittedPaths, hSingleScatterPaths
                if(debug) std::cout << "Got FOMs " << FOMs[i].at(0) << ", " << FOMs[i].at(1) << ", " << FOMs[i].at(2) << std::endl;
                
                bestworstFOMPoints[i] = GetBestFOM(FOMs[i], points[i]);
                if(debug) std::cout << "Got best FOMs" << std::endl;
                
                points[i] = GetThreePoints(bestworstFOMPoints[i].at(0), bestworstFOMPoints[i].at(1), points[i]);
                if(debug) std::cout << "Got new points: " << points[i].at(0) << ", " << points[i].at(1) << ", " << points[i].at(2) << std::endl;
                
                temp_diffs[i] = std::abs(points[i].at(1) - prevBestFOMPoints[i]);
                if(debug) std::cout << "Got difference: " << temp_diffs[i] << std::endl;
                
                prevBestFOMPoints[i] = points[i].at(1);
                if(extraInfo) hFOMxy.at(i)->SetPoint(num_iterations[i], num_iterations[i], bestworstFOMPoints[i].at(2));
                num_iterations[i]++;
            }
            if(debug) std::cout << "Setting fixed point: " << bestworstFOMPoints[i].at(0) << std::endl;
            fixedPoints.at(i) = bestworstFOMPoints[i].at(0);
            if(debug) std::cout << "" << std::endl;
        }

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
        // If any point difference is still larger than its associated tolerance, keep looping.
        loop_condition = false;
        for (int i = 0; i < 6; ++i) {
            if (points_diffs.at(i) > point_tolerances.at(i)) {
                loop_condition = true;
                break;
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
    for (int i = 0; i < 6; ++i) {
        outputFile_txt << std::to_string(fixedPoints.at(i)) + "\n"; //x_a, x_b, x_c, y_a, y_b, y_c
    }
    outputFile_txt.close();

    // now use final result to get region, write to file

    if(verbose) {
        std::cout << "Region found" << std::endl;
        for (int i = 0; i < 6; ++i) {
            std::cout << "    " << point_names[i] << ": " << fixedPoints.at(i) << std::endl;
        }
        std::cout << "There were " << numMainLoopIterations << " iterations of the main loop" << std::endl;
    }

    HistList regionSelectedHists = GetRegionSelectedHists(fixedPoints, hists_lists, fibre, saveroot_txt);

    // get file name from path+filename string
    std::size_t botDirPos = inputFile.find_last_of("/");
    std::string filename = inputFile.substr(botDirPos+1, inputFile.length());
    std::string saveroot = "region_selected_hists_" + signal_param + "_" + filename;
    TFile *rootfile = new TFile(saveroot.c_str(),"RECREATE");
    
    // Write hists to file
    rootfile->cd();
    (TH2F*)regionSelectedHists.Tracking_Hists().at(1)->Write();
    regionSelectedHists.Write();
    
    std::string temp_name;
    if(extraInfo){
        for (int i = 0; i < 6; ++i) {
            temp_name = "hFOM" + point_names[i];
            hFOMxy.at(i)->SetName(temp_name.c_str());
            hFOMxy.at(i)->Write();

            temp_name = "hDiff" + point_names[i];
            hDiffxy.at(i)->SetName(temp_name.c_str());
            hDiffxy.at(i)->Write();

            temp_name = "hPoint" + point_names[i];
            hPointxy.at(i)->SetName(temp_name.c_str());
            hPointxy.at(i)->Write();
        }
        hFOMmainloop->SetName("hFOMmainloop");
        hFOMmainloop->Write();
    }

    // Draw region limits on 2D hist and write to root file
    DrawRegionLims(fixedPoints, hists_lists, fibre);

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
 * @brief Draw region limits as lines on RestHitvsCosTheta histogram.
 * 
 * @param fixedPoints 
 * @param Hists 
 * @param fibre 
 * @return TCanvas* 
 */
void DrawRegionLims(std::vector<double> fixedPoints, HistList hists_lists, std::string fibre) {
    // get points
    double x_a = fixedPoints.at(0);
    double x_b = fixedPoints.at(1);
    double x_c = fixedPoints.at(2);
    double y_a = fixedPoints.at(3);
    double y_b = fixedPoints.at(4);
    double y_c = fixedPoints.at(5);

    double direct_max_time = hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetXaxis()->GetBinCenter(hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetMaximumBin()) + 10;  //hNoEffectPaths
    double direct_min_time = hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetXaxis()->GetBinCenter(hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetMaximumBin()) - 10;  //hNoEffectPaths
    //FIXME: don't hardcode this:
    double direct_cos_alpha;
    if (fibre == "FA089") {  // 10deg off-axis
        direct_cos_alpha = -0.85; 
    } else if (fibre == "FA173" or fibre == "FA150" or fibre == "FA093") {  // 20deg off-axis
        direct_cos_alpha = -0.6;
    } else {  // on-axis
        direct_cos_alpha = -0.9;
    }

    double reflected_max_time = hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetXaxis()->GetBinCenter(hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetMaximumBin()) + 10;  //hNoEffectPaths
    double reflected_min_time = hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetXaxis()->GetBinCenter(hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetMaximumBin()) - 10;  //hNoEffectPaths
    double reflected_cos_alpha = 0.95; //FIXME: don't hardcode this

    // Draw box cuts on resthit vs costheta hist
    TCanvas *c1 = new TCanvas("cuts","cuts");  //Create output canvas to be saved in output file
    TH2F *h = (TH2F*)hists_lists.Tracking_Hists().at(1)->Clone();  //hAllPaths
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

    c1->Write();  
    delete c1;
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

    // Replace the appropriate point in the triangle with each point in points and see if the bin falls in the triangle.
    for (int i = 0; i < 3; ++i){
        Tri[numVar] = points.at(i);
        for(int x=1; x<reEmittedHist->GetNbinsX()+1; x++){ //loop over histogram bins
            double xBinCenter = reEmittedHist->GetXaxis()->GetBinCenter(x);
            for(int y=1; y<reEmittedHist->GetNbinsY()+1; y++){
                double yBinCenter = reEmittedHist->GetYaxis()->GetBinCenter(y);
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
 * @brief Clone all the original histograms 3 times: One for the triangle region,
 * called Region, and one for the direct (called Direct) and reflected (Reflected) beam spots each.
 * The bins falling outside of these regions are then set to zero in each corresponding histogram.
 * Returns vector of resultings histograms.
 * 
 * @param finalPoints final fixed points (x_a_min, x_b_max, x_c_max, y_a_min, y_b_max, y_c_min).
 * @param Hists Original histograms.
 * @param fibre Which fibre fired, to work out regions.
 * @param saveroot_txt Name of txt file to save region lims to.
 * @return std::vector<TH2F*> 
 */
HistList GetRegionSelectedHists(std::vector<double> finalPoints, HistList hists_lists, std::string fibre, std::string saveroot_txt){

    double direct_max_time = hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetXaxis()->GetBinCenter(hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetMaximumBin()) + 10;  //hNoEffectPaths
    double direct_min_time = hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetXaxis()->GetBinCenter(hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetMaximumBin()) - 10;  //hNoEffectPaths
    //FIXME: don't hardcode this:
    double direct_cos_alpha;
    if (fibre == "FA089") {  // 10deg off-axis
        direct_cos_alpha = -0.85; 
    } else if (fibre == "FA173" or fibre == "FA150" or fibre == "FA093") {  // 20deg off-axis
        direct_cos_alpha = -0.6;
    } else {  // on-axis
        direct_cos_alpha = -0.9;
    }

    double reflected_max_time = hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetXaxis()->GetBinCenter(hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetMaximumBin()) + 10;  //hNoEffectPaths
    double reflected_min_time = hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetXaxis()->GetBinCenter(hists_lists.Tracking_Hists().at(5)->ProjectionY()->GetMaximumBin()) - 10;  //hNoEffectPaths
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

    // Create triangle with final points, to check if bin centers are within it in the following loop
    triangle Tri = triangle(finalPoints.at(0), finalPoints.at(1), finalPoints.at(2), finalPoints.at(3),
                            finalPoints.at(4), finalPoints.at(5));

    for(int x=1; x<hists_lists.Tracking_Hists().at(0)->GetNbinsX()+1; x++){ //loop over histogram bins
        double xBinCenter = hists_lists.Tracking_Hists().at(0)->GetXaxis()->GetBinCenter(x);
        for(int y=1; y<hists_lists.Tracking_Hists().at(0)->GetNbinsY()+1; y++){ //loop over histogram bins
            double yBinCenter = hists_lists.Tracking_Hists().at(0)->GetYaxis()->GetBinCenter(y);
            if(!(Tri.check_point_inside_triangle(xBinCenter, yBinCenter))){
                for (int i = 0; i < 15; ++i) {
                    hists_lists.Region_Hists().at(i)->SetBinContent(x,y,0);
                }
            }
            if(xBinCenter > direct_cos_alpha or yBinCenter >= direct_max_time or yBinCenter <= direct_min_time){
                for (int i = 0; i < 15; ++i) {
                    hists_lists.Direct_Hists().at(i)->SetBinContent(x,y,0);
                }
            }
            if(xBinCenter < reflected_cos_alpha or yBinCenter >= reflected_max_time or yBinCenter <= reflected_min_time){
                for (int i = 0; i < 15; ++i) {
                    hists_lists.Reflected_Hists().at(i)->SetBinContent(x,y,0);
                }
            }
        }
    }
    return hists_lists;
}