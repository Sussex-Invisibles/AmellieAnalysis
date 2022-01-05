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

    std::cout << "Step 1" << std::endl;

    TFile *region_hists_file;
    TFile *full_hists_file;

    std::cout << "Step 2" << std::endl;

    try{
        region_hists_file = TFile::Open(region_selected_file.c_str());
        full_hists_file = TFile::Open(full_file.c_str());
    }
    catch(...){
        std::cout << "Could not open input files." << std::endl;
        return 1;
    }

    std::cout << "Step 3" << std::endl;

    if (data_type == "MC") {
        std::cout << "Step 4" << std::endl;

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

        std::cout << "Step 5" << std::endl;

        full_hists_file->GetObject("hPmtResTimeVsCosTheta",fullHist);  // <---------------
        region_hists_file->GetObject("hRegionPmtResTimeVsCosTheta",regionHistAll);  // <--------------------
        region_hists_file->GetObject("hDirectPmtResTimeVsCosTheta",directHistAll);  // <---------------------
        region_hists_file->GetObject("hReflectedPmtResTimeVsCosTheta",reflectedHistAll);
        region_hists_file->GetObject("hRegionReemissionResTimeVsCosTheta",regionHistReemitted);
        region_hists_file->GetObject("hDirectReemissionResTimeVsCosTheta",directHistReemitted);
        region_hists_file->GetObject("hReflectedReemissionResTimeVsCosTheta",reflectedHistReemitted);
        region_hists_file->GetObject("hRegionSingleScatterResTimeVsCosTheta",regionHistScattered);
        region_hists_file->GetObject("hDirectSingleScatterResTimeVsCosTheta",directHistScattered);
        region_hists_file->GetObject("hReflectedSingleScatterResTimeVsCosTheta",reflectedHistScattered);

        std::cout << "Step 6" << std::endl;

        //first stat is signal / signal + background in region

        double sumFullHist = 0;  // <-------------
        double sumFullHistRegion = 0;  // <-------------
        double sumSignalHistRegion = 0;

        double sumFullHistDirect = 0;  // <-------------
        double sumSignalHistDirect = 0;

        double sumFullHistReflected = 0;
        double sumSignalHistReflected = 0;

        std::cout << "Step 7" << std::endl;

        for(int i=0;i<fullHist->GetXaxis()->GetNbins();i++){  // <-------------
            for(int j=0;j<fullHist->GetYaxis()->GetNbins();j++){  // <-------------
                std::cout << "Step 7.1" << std::endl;
                sumFullHist += fullHist->GetBinContent(i,j);  // <-------------
                std::cout << "Step 7.2" << std::endl;
                sumFullHistDirect += directHistAll->GetBinContent(i,j);  // <--------------
                std::cout << "Step 7.3" << std::endl;
                sumFullHistReflected += reflectedHistAll->GetBinContent(i,j);
                std::cout << "Step 7.4" << std::endl;
                sumFullHistRegion += regionHistAll->GetBinContent(i,j);  // <--------------
                std::cout << "Step 7.5" << std::endl;
                if(signal == "reemitted"){
                    std::cout << "Step 7.6" << std::endl;
                    sumSignalHistRegion += regionHistReemitted->GetBinContent(i,j);
                    sumSignalHistDirect += directHistReemitted->GetBinContent(i,j);
                    sumSignalHistReflected += reflectedHistReemitted->GetBinContent(i,j);
                    std::cout << "Step 7.7" << std::endl;
                }
                if(signal == "scattered"){
                    std::cout << "Step 7.8" << std::endl;
                    sumSignalHistRegion += regionHistScattered->GetBinContent(i,j);
                    sumSignalHistDirect += directHistScattered->GetBinContent(i,j);
                    sumSignalHistReflected += reflectedHistScattered->GetBinContent(i,j);
                    std::cout << "Step 7.9" << std::endl;
                }
                if(signal == "attenuated"){
                    std::cout << "Step 7.9.1" << std::endl;
                    sumSignalHistRegion += regionHistScattered->GetBinContent(i,j) + regionHistReemitted->GetBinContent(i,j);
                    sumSignalHistDirect += directHistScattered->GetBinContent(i,j) + directHistReemitted->GetBinContent(i,j);
                    sumSignalHistReflected += reflectedHistScattered->GetBinContent(i,j) + reflectedHistReemitted->GetBinContent(i,j);
                    std::cout << "Step 7.9.2" << std::endl;
                }
            }
        }
        std::cout << "Step 8" << std::endl;

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

        std::cout << "Step 9" << std::endl;

    } else if (data_type == "raw") {
        TH2F *fullHist;  // <---------------
        TH2F *regionHistAll;  // <---------------
        TH2F *directHistAll;  // <---------------
        TH2F *reflectedHistAll;

        full_hists_file->GetObject("hPmtResTimeVsCosTheta",fullHist);  // <---------------
        region_hists_file->GetObject("hRegionPmtResTimeVsCosTheta",regionHistAll);  // <--------------------
        region_hists_file->GetObject("hDirectPmtResTimeVsCosTheta",directHistAll);  // <---------------------
        region_hists_file->GetObject("hReflectedPmtResTimeVsCosTheta",reflectedHistAll);

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

    //~~~ create histograms ~~~//
    HistList hists_lists;
    if (data_type == "MC") {
        hists_lists.Read_File(tracked_file);
    } else if (data_type == "raw") {
        hists_lists.Read_File(tracked_file, {"hPmtResTimeVsCosTheta"});
    } else {
        std::cout << "Wrong data type. Should be MC or raw" << std::endl;
        throw;
    }
        
    //~~~ Apply region cuts by looping through bins ~~~//

    //direct beam spot, max +/- 10 ns, -0.9 cos theta (straight through)
    double min_angle_direct_beam_spot = -0.9;
    //reflected beam spot, max +/- 10 ns, 0.95 cos theta (straight through)
    double min_angle_reflected_beam_spot = 0.95;

    int nBinsX = hists_lists.Tracking_Hists().at(0)->GetXaxis()->GetNbins();
    int nBinsY = hists_lists.Tracking_Hists().at(0)->GetYaxis()->GetNbins();
    double xBinCenter;
    double yBinCenter;
    bool direct_Xcut;
    bool reflected_Xcut;
    for (int x=0; x<nBinsX+1; x++) {
        xBinCenter = hists_lists.Tracking_Hists().at(0)->GetXaxis()->GetBinCenter(x);
        if (xBinCenter >= min_angle_direct_beam_spot) {direct_Xcut = true;} else {direct_Xcut = false;}
        if (xBinCenter <= min_angle_reflected_beam_spot) {reflected_Xcut = true;} else {reflected_Xcut = false;}
        for (int y=0; y<nBinsY+1; y++) {
            yBinCenter = hists_lists.Tracking_Hists().at(0)->GetYaxis()->GetBinCenter(y);
            // Check if bin is within triangle region
            if (!(Tri.check_point_inside_triangle(xBinCenter, yBinCenter))) {
                for (int i = 0; i < 15; ++i) {
                    hists_lists.Region_Hists().at(i)->SetBinContent(x, y, 0);
                }
            }
            // Check if bin is in direct beam spot
            if (direct_Xcut) {
                if (yBinCenter <= min_time_direct_beam_spot or yBinCenter >= max_time_direct_beam_spot) {
                    for (int i = 0; i < 15; ++i) {
                        hists_lists.Direct_Hists().at(i)->SetBinContent(x, y, 0);
                    }
                }
            }
            // Check if bin is in direct beam spot
            if (reflected_Xcut) {
                if (yBinCenter <= min_time_reflected_beam_spot or yBinCenter >= max_time_reflected_beam_spot) {
                    for (int i = 0; i < 15; ++i) {
                        hists_lists.Reflected_Hists().at(i)->SetBinContent(x, y, 0);
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
    hists_lists.Write();
    rootfile->Write();
    rootfile->Close();

    return 0;
}