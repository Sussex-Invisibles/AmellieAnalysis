//Compile: g++ -g -std=c++1y -o generateStats.exe generateStats.cpp `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TVectorD.h>
#include "AMELLIE_utils.hpp"

/*
Either: 
    - supply a region selected file and generate stats (FOMs)
    - supply region co-ordinates and generate a region selected file
*/

int generate_stats(std::string region_selected_file, std::string full_file, std::string signal, std::string data_type);
int make_region_cut(std::string tracked_file, triangle Tri, double min_time_direct_beam_spot, double max_time_direct_beam_spot,
                    double min_time_reflected_beam_spot, double max_time_reflected_beam_spot, std::string data_type);

int main(int argv, char** argc){
    std::string choice = argc[1];
    if(choice == "generate_stats"){
        std::cout << "Generating stats" << std::endl;
        std::string region_selected_file = argc[2];
        std::string full_file = argc[3];
        std::string signal = argc[4];
        std::string data_type = argc[5];  // MC or raw
        int status = generate_stats(region_selected_file, full_file, signal, data_type);
        return status;
    }
    else if(choice == "apply_region"){
        std::cout << "Applying region cuts" << std::endl;
        std::string tracked_file = argc[2];
        triangle Tri = triangle(std::stod(argc[3]), std::stod(argc[4]), std::stod(argc[5]), std::stod(argc[6]),
                                std::stod(argc[7]), std::stod(argc[8]));
        double min_direct = std::stod(argc[9]);
        double max_direct = std::stod(argc[10]);
        double min_reflected = std::stod(argc[11]);
        double max_reflected = std::stod(argc[12]);
        std::string data_type = argc[13];  // MC or raw
        int status = make_region_cut(tracked_file, Tri, min_direct, max_direct, min_reflected, max_reflected, data_type);
        return status;
    }
    else{
        std::cout << "Unknown input, use either generate_stats or apply_region" << std::endl;
        return 1;
    }
}


/**
 * @brief Create stats on the different regions (triangle and direct and reflected beam spots):
 * sumSignalHistRegion / sqrt(sumFullHistRegion),
 * sumFullHist,
 * sumFullHistRegion,
 * sumSignalHistRegion,
 * sumFullHistDirect,
 * sumSignalHistDirect,
 * sumFullHistReflected,
 * sumSignalHistReflected.
 * 
 * @param region_selected_file root file with region histograms.
 * @param full_file original file with overall histograms.
 * @param signal = reemitted, scattered or attenuated(=reemitted+scattered). Signal that was optimised for.
 * @return int 
 */
int generate_stats(std::string region_selected_file, std::string full_file, std::string signal, std::string data_type){
    // read in files

    if(signal != "reemitted" and signal != "scattered" and signal != "attenuated"){ // attenuated = scattered + reemitted
        std::cout << "Invalid signal! Use: reemitted, scattered, attenuated";
        return 1;
    }

    TFile *region_hists_file;
    TFile *full_hists_file;

    try{
        region_hists_file = TFile::Open(region_selected_file.c_str());
        full_hists_file = TFile::Open(full_file.c_str());
    }
    catch(...){
        std::cout << "Could not open input files." << std::endl;
        return 1;
    }

    if (data_type == "MC") {
        TH2F *fullHist;  // <---------------
        TH2F *regionHistReemitted;
        TH2F *directHistReemitted;
        TH2F *reflectedHistReemitted;
        TH2F *regionHistScattered;
        TH2F *directHistScattered;
        TH2F *reflectedHistScattered;
        TH2F *regionHistAll;  // <---------------
        TH2F *directHistAll;  // <---------------
        TH2F *reflectedHistAll;

        full_hists_file->GetObject("hPmtResTimeVsCosTheta",fullHist);  // <---------------
        region_hists_file->GetObject("hRegionSelectedReEmittedPaths",regionHistReemitted);
        region_hists_file->GetObject("hDirectCutReEmittedPaths",directHistReemitted);
        region_hists_file->GetObject("hReflectedCutReEmittedPaths",reflectedHistReemitted);
        region_hists_file->GetObject("hRegionCutSingleScatterPaths",regionHistScattered);
        region_hists_file->GetObject("hDirectCutSingleScatterPaths",directHistScattered);
        region_hists_file->GetObject("hReflectedCutSingleScatterPaths",reflectedHistScattered);
        region_hists_file->GetObject("hRegionSelectedAllPaths",regionHistAll);  // <--------------------
        region_hists_file->GetObject("hDirectCutAllPaths",directHistAll);  // <---------------------
        region_hists_file->GetObject("hReflectedCutAllPaths",reflectedHistAll);

        //first stat is signal / signal + background in region

        double sumFullHist = 0;  // <-------------
        double sumFullHistRegion = 0;  // <-------------
        double sumSignalHistRegion = 0;

        double sumFullHistDirect = 0;  // <-------------
        double sumSignalHistDirect = 0;

        double sumFullHistReflected = 0;
        double sumSignalHistReflected = 0;

        for(int i=0;i<fullHist->GetXaxis()->GetNbins();i++){  // <-------------
            for(int j=0;j<fullHist->GetYaxis()->GetNbins();j++){  // <-------------
                sumFullHist += fullHist->GetBinContent(i,j);  // <-------------
                sumFullHistDirect += directHistAll->GetBinContent(i,j);  // <--------------
                sumFullHistReflected += reflectedHistAll->GetBinContent(i,j);
                sumFullHistRegion += regionHistAll->GetBinContent(i,j);  // <--------------
                if(signal == "reemitted"){
                    sumSignalHistRegion += regionHistReemitted->GetBinContent(i,j);
                    sumSignalHistDirect += directHistReemitted->GetBinContent(i,j);
                    sumSignalHistReflected += reflectedHistReemitted->GetBinContent(i,j);
                }
                if(signal == "scattered"){
                    sumSignalHistRegion += regionHistScattered->GetBinContent(i,j);
                    sumSignalHistDirect += directHistScattered->GetBinContent(i,j);
                    sumSignalHistReflected += reflectedHistScattered->GetBinContent(i,j);
                }
                if(signal == "attenuated"){
                    sumSignalHistRegion += regionHistScattered->GetBinContent(i,j) + regionHistReemitted->GetBinContent(i,j);
                    sumSignalHistDirect += directHistScattered->GetBinContent(i,j) + directHistReemitted->GetBinContent(i,j);
                    sumSignalHistReflected += reflectedHistScattered->GetBinContent(i,j) + reflectedHistReemitted->GetBinContent(i,j);
                }
            }
        }

        // open txt file and write values, FOMs can be made and plotted later!

        std::ofstream outputFile;
        // get file name from path+filename string
        std::size_t region_botDirPos = region_selected_file.find_last_of("/");
        std::string region_filename = region_selected_file.substr(region_botDirPos+1, region_selected_file.length());
        std::string saveroot = "stats_for_" + region_filename.substr(0, region_filename.length()-5) + ".txt";
        outputFile.open(saveroot.c_str());
        outputFile << signal + "\n";
        outputFile << std::to_string(sumSignalHistRegion / sqrt(sumFullHistRegion)) + "\n";
        outputFile << std::to_string(sumFullHist) + "\n";  // <-------------
        outputFile << std::to_string(sumFullHistRegion) + "\n";  // <-------------
        outputFile << std::to_string(sumSignalHistRegion) + "\n";
        outputFile << std::to_string(sumFullHistDirect) + "\n";  // <-------------
        outputFile << std::to_string(sumSignalHistDirect) + "\n";
        outputFile << std::to_string(sumFullHistReflected) + "\n";
        outputFile << std::to_string(sumSignalHistReflected) + "\n";
        outputFile.close();

    } else if (data_type == "raw") {
        TH2F *fullHist;  // <---------------
        TH2F *regionHistAll;  // <---------------
        TH2F *directHistAll;  // <---------------
        TH2F *reflectedHistAll;

        full_hists_file->GetObject("hPmtResTimeVsCosTheta",fullHist);  // <---------------
        region_hists_file->GetObject("hRegionSelectedAllPaths",regionHistAll);  // <--------------------
        region_hists_file->GetObject("hDirectCutAllPaths",directHistAll);  // <---------------------
        region_hists_file->GetObject("hReflectedCutAllPaths",reflectedHistAll);

        //first stat is signal / signal + background in region

        double sumFullHist = 0;
        double sumFullHistRegion = 0;
        double sumFullHistDirect = 0;
        double sumFullHistReflected = 0;

        for(int i=0;i<fullHist->GetXaxis()->GetNbins();i++){  // <-------------
            for(int j=0;j<fullHist->GetYaxis()->GetNbins();j++){  // <-------------
                sumFullHist += fullHist->GetBinContent(i,j);  // <-------------
                sumFullHistDirect += directHistAll->GetBinContent(i,j);  // <--------------
                sumFullHistReflected += reflectedHistAll->GetBinContent(i,j);
                sumFullHistRegion += regionHistAll->GetBinContent(i,j);  // <--------------
            }
        }

        // open txt file and write values, FOMs can be made and plotted later!

        std::ofstream outputFile;
        // get file name from path+filename string
        std::size_t region_botDirPos = region_selected_file.find_last_of("/");
        std::string region_filename = region_selected_file.substr(region_botDirPos+1, region_selected_file.length());
        std::string saveroot = "stats_for_" + region_filename.substr(0, region_filename.length()-5) + ".txt";
        outputFile.open(saveroot.c_str());
        outputFile << signal + "\n";
        outputFile << std::to_string(sumFullHist) + "\n";  // <-------------
        outputFile << std::to_string(sumFullHistRegion) + "\n";  // <-------------
        outputFile << std::to_string(sumFullHistDirect) + "\n";  // <-------------
        outputFile << std::to_string(sumFullHistReflected) + "\n";
        outputFile.close();

    } else {
        std::cout << "Wrong data type. Should be MC or raw" << std::endl;
        throw;
    }

    return 0;
}

/**
 * @brief Make root file with region histograms, just as getRegions.cpp does.
 * Points: a - bottom-left, b - bottom-right, c - top right. (y-axis = residual
 * time, x-axis = cos(angle from AV centre to fibre direction)).
 * x-limits of direct/reflected beam spot limits are hardcoded in the function.
 * 
 * @param tracked_file Original root file with total histograms.
 * @param Tri triangle object with coordinates of each vertex saved.
 * @param min_time_direct_beam_spot Direct beam spot limit.
 * @param max_time_direct_beam_spot Direct beam spot limit.
 * @param min_time_reflected_beam_spot Reflected beam spot limit.
 * @param max_time_reflected_beam_spot Reflected beam spot limit.
 * @return int 
 */
int make_region_cut(std::string tracked_file, triangle Tri,
                    double min_time_direct_beam_spot, double max_time_direct_beam_spot, double min_time_reflected_beam_spot,
                    double max_time_reflected_beam_spot, std::string data_type){

    TFile *tracked_hists_file;

    try{
        tracked_hists_file = TFile::Open(tracked_file.c_str());
    }
    catch(...){
        std::cout << "Could not open input files." << std::endl;
        return 1;
    }

    if (data_type == "MC") {
        //~~~ create histograms and graphs ~~~//

        // Tracking root file hists
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
        std::vector<const char*> Hist_names = {"hReemissionResTimeVsCosTheta", "hPmtResTimeVsCosTheta", "hNoiseResTimeVsCosTheta", 
                                    "hSingleScatterResTimeVsCosTheta", "hOtherEffectResTimeVsCosTheta", "hNoEffectResTimeVsCosTheta", 
                                    "hNearReflectResTimeVsCosTheta", "hRopesResTimeVsCosTheta", "hPMTReflectionResTimeVsCosTheta", 
                                    "hExtWaterScatterResTimeVsCosTheta", "hInnerAvReflectionResTimeVsCosTheta", "hAVPipesResTimeVsCosTheta", 
                                    "hAcrylicScatterResTimeVsCosTheta", "OtherScatterResTimeVsCosTheta"};


        // New hists to clone to
        TH2F *hRegionSelectedReEmittedPaths; TH2F *hRegionSelectedAllPaths;
        TH2F *hRegionSelectedNoisePaths; TH2F *hRegionSelectedSingleScatterPaths;
        TH2F *hRegionSelectedOtherPaths; TH2F *hRegionSelectedNoEffectPaths;
        TH2F *hRegionSelectedNearReflectPaths; TH2F *hRegionSelectedRopesPaths;
        TH2F *hRegionSelectedPMTReflectionPaths; TH2F *hRegionSelectedExtWaterScatterPaths;
        TH2F *hRegionSelectedInnerAvReflectPaths; TH2F *hRegionSelectedMultipleEffectPaths;
        TH2F *hRegionSelectedAVPipesPaths; TH2F *hRegionSelectedAcrylicPaths;
        TH2F *hRegionSelectedOtherScatterPaths;
        std::vector<TH2F*> RegionHists;
        RegionHists.push_back(hRegionSelectedReEmittedPaths); RegionHists.push_back(hRegionSelectedAllPaths);
        RegionHists.push_back(hRegionSelectedNoisePaths); RegionHists.push_back(hRegionSelectedSingleScatterPaths);
        RegionHists.push_back(hRegionSelectedOtherPaths); RegionHists.push_back(hRegionSelectedNoEffectPaths);
        RegionHists.push_back(hRegionSelectedNearReflectPaths); RegionHists.push_back(hRegionSelectedRopesPaths);
        RegionHists.push_back(hRegionSelectedPMTReflectionPaths); RegionHists.push_back(hRegionSelectedExtWaterScatterPaths);
        RegionHists.push_back(hRegionSelectedInnerAvReflectPaths); RegionHists.push_back(hRegionSelectedMultipleEffectPaths);
        RegionHists.push_back(hRegionSelectedAVPipesPaths); RegionHists.push_back(hRegionSelectedAcrylicPaths);
        RegionHists.push_back(hRegionSelectedOtherScatterPaths);
        std::vector<const char*> RegionHist_names = {"hRegionSelectedReEmittedPaths", "hRegionSelectedAllPaths", "hRegionCutNoisePaths",
                                            "hRegionCutSingleScatterPaths", "hRegionCutOtherPaths", "hRegionCutNoEffectPaths"
                                            "hRegionCutNearReflectPaths", "hRegionCutRopesPaths", "hRegionCutPMTReflectionPaths"
                                            "hRegionCutExtWaterScatterPaths", "hRegionCutInnerAvReflectPaths", "hRegionCutMultipleEffectPaths"
                                            "hRegionCutAVPipesPaths", "hRegionCutAcrylicPaths", "hRegionCutOtherScatterPaths"};

        TH2F *hDirectCutReEmittedPaths; TH2F *hDirectCutAllPaths;
        TH2F *hDirectCutNoisePaths; TH2F *hDirectCutSingleScatterPaths;
        TH2F *hDirectCutOtherPaths; TH2F *hDirectCutNoEffectPaths;
        TH2F *hDirectCutNearReflectPaths; TH2F *hDirectCutRopesPaths;
        TH2F *hDirectCutPMTReflectionPaths; TH2F *hDirectCutExtWaterScatterPaths;
        TH2F *hDirectCutInnerAvReflectPaths; TH2F *hDirectCutMultipleEffectPaths;
        TH2F *hDirectCutAVPipesPaths; TH2F *hDirectCutAcrylicPaths;
        TH2F *hDirectCutOtherScatterPaths;
        std::vector<TH2F*> DirectHists;
        DirectHists.push_back(hDirectCutReEmittedPaths); DirectHists.push_back(hDirectCutAllPaths);
        DirectHists.push_back(hDirectCutNoisePaths); DirectHists.push_back(hDirectCutSingleScatterPaths);
        DirectHists.push_back(hDirectCutOtherPaths); DirectHists.push_back(hDirectCutNoEffectPaths);
        DirectHists.push_back(hDirectCutNearReflectPaths); DirectHists.push_back(hDirectCutRopesPaths);
        DirectHists.push_back(hDirectCutPMTReflectionPaths); DirectHists.push_back(hDirectCutExtWaterScatterPaths);
        DirectHists.push_back(hDirectCutInnerAvReflectPaths); DirectHists.push_back(hDirectCutMultipleEffectPaths);
        DirectHists.push_back(hDirectCutAVPipesPaths); DirectHists.push_back(hDirectCutAcrylicPaths);
        DirectHists.push_back(hDirectCutOtherScatterPaths);
        std::vector<const char*> DirectHist_names = {"hDirectCutReEmittedPaths", "hDirectCutAllPaths", "hDirectCutNoisePaths", 
                                            "hDirectCutSingleScatterPaths", "hDirectCutOtherPaths", "hDirectCutNoEffectPaths", 
                                            "hDirectCutNearReflectPaths", "hDirectCutRopesPaths", "hDirectCutPMTReflectionPaths", 
                                            "hDirectCutExtWaterScatterPaths", "hDirectCutInnerAvReflectPaths", "hDirectCutAVPipesPaths", 
                                            "hDirectCutAcrylicPaths", "hDirectCutOtherScatterPaths"};

        TH2F *hReflectedCutReEmittedPaths; TH2F *hReflectedCutAllPaths;
        TH2F *hReflectedCutNoisePaths; TH2F *hReflectedCutSingleScatterPaths;
        TH2F *hReflectedCutOtherPaths; TH2F *hReflectedCutNoEffectPaths;
        TH2F *hReflectedCutNearReflectPaths; TH2F *hReflectedCutRopesPaths;
        TH2F *hReflectedCutPMTReflectionPaths; TH2F *hReflectedCutExtWaterScatterPaths;
        TH2F *hReflectedCutInnerAvReflectPaths; TH2F *hReflectedCutMultipleEffectPaths;
        TH2F *hReflectedCutAVPipesPaths; TH2F *hReflectedCutAcrylicPaths;
        TH2F *hReflectedCutOtherScatterPaths;
        std::vector<TH2F*> ReflectedHists;
        ReflectedHists.push_back(hReflectedCutReEmittedPaths); ReflectedHists.push_back(hReflectedCutAllPaths);
        ReflectedHists.push_back(hReflectedCutNoisePaths); ReflectedHists.push_back(hReflectedCutSingleScatterPaths);
        ReflectedHists.push_back(hReflectedCutOtherPaths); ReflectedHists.push_back(hReflectedCutNoEffectPaths);
        ReflectedHists.push_back(hReflectedCutNearReflectPaths); ReflectedHists.push_back(hReflectedCutRopesPaths);
        ReflectedHists.push_back(hReflectedCutPMTReflectionPaths); ReflectedHists.push_back(hReflectedCutExtWaterScatterPaths);
        ReflectedHists.push_back(hReflectedCutInnerAvReflectPaths); ReflectedHists.push_back(hReflectedCutMultipleEffectPaths);
        ReflectedHists.push_back(hReflectedCutAVPipesPaths); ReflectedHists.push_back(hReflectedCutAcrylicPaths);
        ReflectedHists.push_back(hReflectedCutOtherScatterPaths);
        std::vector<const char*> ReflectedHist_names = {"hReflectedCutReEmittedPaths", "hReflectedCutAllPaths", "hReflectedCutNoisePaths", 
                                            "hReflectedCutSingleScatterPaths", "hReflectedCutOtherPaths", "hReflectedCutNoEffectPaths", 
                                            "hReflectedCutNearReflectPaths", "hReflectedCutRopesPaths", "hReflectedCutPMTReflectionPaths", 
                                            "hReflectedCutExtWaterScatterPaths", "hReflectedCutInnerAvReflectPaths", "hReflectedCutMultipleEffectPaths", 
                                            "hReflectedCutAVPipesPaths", "hReflectedCutAcrylicPaths", "hReflectedCutOtherScatterPaths"};

        // Assign hists
        for (int i = 0; i < 15; ++i) {
            tracked_hists_file->GetObject(Hist_names.at(i), Hists.at(i));

            RegionHists.at(i) = (TH2F*)Hists.at(i)->Clone();
            RegionHists.at(i)->SetName(RegionHist_names.at(i));

            DirectHists.at(i) = (TH2F*)Hists.at(i)->Clone();
            DirectHists.at(i)->SetName(DirectHist_names.at(i));

            ReflectedHists.at(i) = (TH2F*)Hists.at(i)->Clone();
            ReflectedHists.at(i)->SetName(ReflectedHist_names.at(i));
        }

        //~~~ Apply region cuts by looping through bins ~~~//

        //direct beam spot, max +/- 10 ns, -0.9 cos theta (straight through)
        double min_angle_direct_beam_spot = -0.9;
        //reflected beam spot, max +/- 10 ns, 0.95 cos theta (straight through)
        double min_angle_reflected_beam_spot = 0.95;

        int nBinsX = hAllPaths->GetXaxis()->GetNbins();
        int nBinsY = hAllPaths->GetYaxis()->GetNbins();
        double xBinCenter;
        double yBinCenter;
        bool direct_Xcut;
        bool reflected_Xcut;
        for (int x=0; x<nBinsX+1; x++) {
            xBinCenter = hAllPaths->GetXaxis()->GetBinCenter(x);
            if (xBinCenter >= min_angle_direct_beam_spot) {direct_Xcut = true;} else {direct_Xcut = false;}
            if (xBinCenter >= min_angle_reflected_beam_spot) {reflected_Xcut = true;} else {reflected_Xcut = false;}
            for (int y=0; y<nBinsY+1; y++) {
                yBinCenter = hAllPaths->GetYaxis()->GetBinCenter(y);

                // Check if bin is within triangle region
                if (!(Tri.check_point_inside_triangle(xBinCenter, yBinCenter))) {
                    for (int i = 0; i < 15; ++i) {
                        RegionHists.at(i)->SetBinContent(x, y, 0);
                    }
                }

                // Check if bin is in direct beam spot
                if (direct_Xcut) {
                    if (yBinCenter <= min_time_direct_beam_spot or yBinCenter >= max_time_direct_beam_spot) {
                        for (int i = 0; i < 15; ++i) {
                            DirectHists.at(i)->SetBinContent(x, y, 0);
                        }
                    }
                }

                // Check if bin is in direct beam spot
                if (reflected_Xcut) {
                    if (yBinCenter <= min_time_reflected_beam_spot or yBinCenter >= max_time_reflected_beam_spot) {
                        for (int i = 0; i < 15; ++i) {
                            DirectHists.at(i)->SetBinContent(x, y, 0);
                        }
                    }
                }
            }
        }

        //write to file
        // get file name from path+filename string
        std::size_t tracked_botDirPos = tracked_file.find_last_of("/");
        std::string tracked_filename = tracked_file.substr(tracked_botDirPos+1, tracked_file.length());
        std::string saveroot = "region_selected_hists_x_a" + std::to_string(Tri.X_a())  + "_x_b_" + std::to_string(Tri.X_b()) + "_x_c_"
                                + std::to_string(Tri.X_c()) + "_y_a_" + std::to_string(Tri.Y_a()) + "_y_b_" + std::to_string(Tri.Y_b())
                                + "_y_c_" + std::to_string(Tri.Y_c()) + "_" + tracked_filename;
        TFile *rootfile = new TFile(saveroot.c_str(),"RECREATE");

        rootfile->cd();
        for (int i = 0; i < 15; ++i) {
            RegionHists.at(i)->Write();
            DirectHists.at(i)->Write();
            ReflectedHists.at(i)->Write();
        }
        rootfile->Write();
        rootfile->Close();

    } else if (data_type == "raw") {
        //~~~ create histograms and graphs ~~~//

        // Tracking root file hists
        TH2F *hAllPaths;
        tracked_hists_file->GetObject("hPmtResTimeVsCosTheta",hAllPaths);

        // New hists to clone to
        TH2F *hRegionSelectedAllPaths = (TH2F*)hAllPaths->Clone();
        hRegionSelectedAllPaths->SetName("hRegionSelectedAllPaths");
        TH2F *hDirectCutAllPaths = (TH2F*)hAllPaths->Clone();
        hDirectCutAllPaths->SetName("hDirectCutAllPaths");
        TH2F *hReflectedCutAllPaths = (TH2F*)hAllPaths->Clone();
        hReflectedCutAllPaths->SetName("hReflectedCutAllPaths");

        //~~~ Apply region cuts by looping through bins ~~~//

        //direct beam spot, max +/- 10 ns, -0.9 cos theta (straight through)
        double min_angle_direct_beam_spot = -0.9;
        //reflected beam spot, max +/- 10 ns, 0.95 cos theta (straight through)
        double min_angle_reflected_beam_spot = 0.95;

        int nBinsX = hAllPaths->GetXaxis()->GetNbins();
        int nBinsY = hAllPaths->GetYaxis()->GetNbins();
        double xBinCenter;
        double yBinCenter;
        bool direct_Xcut;
        bool reflected_Xcut;
        for (int x=0; x<nBinsX+1; x++) {
            xBinCenter = hAllPaths->GetXaxis()->GetBinCenter(x);
            if (xBinCenter >= min_angle_direct_beam_spot) {direct_Xcut = true;} else {direct_Xcut = false;}
            if (xBinCenter >= min_angle_reflected_beam_spot) {reflected_Xcut = true;} else {reflected_Xcut = false;}
            for (int y=0; y<nBinsY+1; y++) {
                yBinCenter = hAllPaths->GetYaxis()->GetBinCenter(y);

                // Check if bin is within triangle region
                if (!(Tri.check_point_inside_triangle(xBinCenter, yBinCenter))) {
                    hRegionSelectedAllPaths->SetBinContent(x, y, 0);
                }
                // Check if bin is in direct beam spot
                if (direct_Xcut) {
                    if (yBinCenter <= min_time_direct_beam_spot or yBinCenter >= max_time_direct_beam_spot) {
                        hDirectCutAllPaths->SetBinContent(x, y, 0);
                    }
                }
                // Check if bin is in direct beam spot
                if (reflected_Xcut) {
                    if (yBinCenter <= min_time_reflected_beam_spot or yBinCenter >= max_time_reflected_beam_spot) {
                        hReflectedCutAllPaths->SetBinContent(x, y, 0);
                    }
                }
            }
        }

        //write to file
        // get file name from path+filename string
        std::size_t tracked_botDirPos = tracked_file.find_last_of("/");
        std::string tracked_filename = tracked_file.substr(tracked_botDirPos+1, tracked_file.length());
        std::string saveroot = "region_selected_hists_x_a_" + std::to_string(Tri.X_a())  + "_x_b_" + std::to_string(Tri.X_b())
                                + "_x_c_" + std::to_string(Tri.X_c()) + "_y_a_" + std::to_string(Tri.Y_a()) + "_y_b_" + std::to_string(Tri.Y_b())
                                + "_y_c_" + std::to_string(Tri.Y_c()) + "_" + tracked_filename;
        TFile *rootfile = new TFile(saveroot.c_str(),"RECREATE");

        rootfile->cd();
        hRegionSelectedAllPaths->Write();
        hDirectCutAllPaths->Write();
        hReflectedCutAllPaths->Write();
        rootfile->Write();
        rootfile->Close();

    } else {
        std::cout << "Wrong data type. Should be MC or raw" << std::endl;
        throw;
    }

    return 0;
}