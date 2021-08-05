//Compile: g++ -g -std=c++1y -o generateStats.exe generateStats.cpp `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TVectorD.h>

/*
Either: 
    - supply a region selected file and generate stats (FOMs)
    - supply region co-ordinates and generate a region selected file
*/

int generate_stats(std::string region_selected_file, std::string full_file, std::string signal);
int make_region_cut(std::string tracked_file, double x_a, double x_b, double x_c, double y_a, double y_b, double y_c, double min_time_direct_beam_spot, double max_time_direct_beam_spot, double min_time_reflected_beam_spot, double max_time_reflected_beam_spot);

int main(int argv, char** argc){
    std::string choice = argc[1];
    if(choice == "generate_stats"){
        std::cout << "Generating stats" << std::endl;
        std::string region_selected_file = argc[2];
        std::string full_file = argc[3];
        std::string signal = argc[4];
        int status = generate_stats(region_selected_file, full_file, signal);
        return status;
    }
    else if(choice == "apply_region"){
        std::cout << "Applying region cuts" << std::endl;
        std::string tracked_file = argc[2];
        double x_a = std::stod(argc[3]);
        double x_b = std::stod(argc[4]);
        double x_c = std::stod(argc[5]);
        double y_a = std::stod(argc[6]);
        double y_b = std::stod(argc[7]);
        double y_c = std::stod(argc[8]);
        double min_direct = std::stod(argc[9]);
        double max_direct = std::stod(argc[10]);
        double min_reflected = std::stod(argc[11]);
        double max_reflected = std::stod(argc[12]);
        int status = make_region_cut(tracked_file, x_a, x_b, x_c, y_a, y_b, y_c, min_direct, max_direct, min_reflected, max_reflected);
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
int generate_stats(std::string region_selected_file, std::string full_file, std::string signal){
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

    TH2F *fullHist;
    TH2F *regionHistReemitted;
    TH2F *directHistReemitted;
    TH2F *reflectedHistReemitted;
    TH2F *regionHistScattered;
    TH2F *directHistScattered;
    TH2F *reflectedHistScattered;
    TH2F *regionHistAll;
    TH2F *directHistAll;
    TH2F *reflectedHistAll;

    full_hists_file->GetObject("hPmtResTimeVsCosTheta",fullHist);
    region_hists_file->GetObject("hRegionSelectedReEmittedPaths",regionHistReemitted);
    region_hists_file->GetObject("hDirectCutReEmittedPaths",directHistReemitted);
    region_hists_file->GetObject("hReflectedCutReEmittedPaths",reflectedHistReemitted);
    region_hists_file->GetObject("hRegionCutSingleScatterPaths",regionHistScattered);
    region_hists_file->GetObject("hDirectCutSingleScatterPaths",directHistScattered);
    region_hists_file->GetObject("hReflectedCutSingleScatterPaths",reflectedHistScattered);
    region_hists_file->GetObject("hRegionSelectedAllPaths",regionHistAll);
    region_hists_file->GetObject("hDirectCutAllPaths",directHistAll);
    region_hists_file->GetObject("hReflectedCutAllPaths",reflectedHistAll);

    //first stat is signal / signal + background in region

    double sumFullHist = 0;
    double sumFullHistRegion = 0;
    double sumSignalHistRegion = 0;

    double sumFullHistDirect = 0;
    double sumSignalHistDirect = 0;

    double sumFullHistReflected = 0;
    double sumSignalHistReflected = 0;

    for(int i=0;i<fullHist->GetXaxis()->GetNbins();i++){
        for(int j=0;j<fullHist->GetYaxis()->GetNbins();j++){
            sumFullHist = sumFullHist + fullHist->GetBinContent(i,j);
            sumFullHistDirect = sumFullHistDirect + directHistAll->GetBinContent(i,j);
            sumFullHistReflected = sumFullHistReflected + reflectedHistAll->GetBinContent(i,j);
            sumFullHistRegion = sumFullHistRegion + regionHistAll->GetBinContent(i,j);
            if(signal == "reemitted"){
                sumSignalHistRegion = sumSignalHistRegion + regionHistReemitted->GetBinContent(i,j);
                sumSignalHistDirect = sumSignalHistDirect + directHistReemitted->GetBinContent(i,j);
                sumSignalHistReflected = sumSignalHistReflected + reflectedHistReemitted->GetBinContent(i,j);
            }
            if(signal == "scattered"){
                sumSignalHistRegion = sumSignalHistRegion + regionHistScattered->GetBinContent(i,j);
                sumSignalHistDirect = sumSignalHistDirect + directHistScattered->GetBinContent(i,j);
                sumSignalHistReflected = sumSignalHistReflected + reflectedHistScattered->GetBinContent(i,j);
            }
            if(signal == "attenuated"){
                sumSignalHistRegion = sumSignalHistRegion + regionHistScattered->GetBinContent(i,j) + regionHistReemitted->GetBinContent(i,j);
                sumSignalHistDirect = sumSignalHistDirect + directHistScattered->GetBinContent(i,j) + directHistReemitted->GetBinContent(i,j);
                sumSignalHistReflected = sumSignalHistReflected + reflectedHistScattered->GetBinContent(i,j) + reflectedHistReemitted->GetBinContent(i,j);
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
    outputFile << std::to_string(sumFullHist) + "\n";
    outputFile << std::to_string(sumFullHistRegion) + "\n";
    outputFile << std::to_string(sumSignalHistRegion) + "\n";
    outputFile << std::to_string(sumFullHistDirect) + "\n";
    outputFile << std::to_string(sumSignalHistDirect) + "\n";
    outputFile << std::to_string(sumFullHistReflected) + "\n";
    outputFile << std::to_string(sumSignalHistReflected) + "\n";
    outputFile.close();

    return 0;
}

/**
 * @brief Make root file with region histograms, just as getRegions.cpp does.
 * Points: a - bottom-left, b - bottom-right, c - top right. (y-axis = residual
 * time, x-axis = cos(angle from AV centre to fibre direction)).
 * x-limits of direct/reflected beam spot limits are hardcoded in the function.
 * 
 * @param tracked_file Original root file with total histograms.
 * @param x_a Triangle point.
 * @param x_b Triangle point.
 * @param x_c Triangle point.
 * @param y_a Triangle point.
 * @param y_b Triangle point.
 * @param y_c Triangle point.
 * @param min_time_direct_beam_spot Direct beam spot limit.
 * @param max_time_direct_beam_spot Direct beam spot limit.
 * @param min_time_reflected_beam_spot Reflected beam spot limit.
 * @param max_time_reflected_beam_spot Reflected beam spot limit.
 * @return int 
 */
int make_region_cut(std::string tracked_file, double x_a, double x_b, double x_c, double y_a, double y_b, double y_c, double min_time_direct_beam_spot, double max_time_direct_beam_spot, double min_time_reflected_beam_spot, double max_time_reflected_beam_spot){

    TFile *tracked_hists_file;

    try{
        tracked_hists_file = TFile::Open(tracked_file.c_str());
    }
    catch(...){
        std::cout << "Could not open input files." << std::endl;
        return 1;
    }

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

    tracked_hists_file->GetObject("hReemissionResTimeVsCosTheta",hReEmittedPaths);
    tracked_hists_file->GetObject("hPmtResTimeVsCosTheta",hAllPaths);

    tracked_hists_file->GetObject("hNoiseResTimeVsCosTheta",hNoisePaths);
    tracked_hists_file->GetObject("hSingleScatterResTimeVsCosTheta",hSingleScatterPaths);
    tracked_hists_file->GetObject("hOtherEffectResTimeVsCosTheta",hOtherPaths);
    tracked_hists_file->GetObject("hNoEffectResTimeVsCosTheta",hNoEffectPaths);
    tracked_hists_file->GetObject("hNearReflectResTimeVsCosTheta",hNearReflectPaths);
    tracked_hists_file->GetObject("hRopesResTimeVsCosTheta",hRopesPaths);
    tracked_hists_file->GetObject("hPMTReflectionResTimeVsCosTheta",hPMTReflectionPaths);
    tracked_hists_file->GetObject("hExtWaterScatterResTimeVsCosTheta",hExtWaterScatterPaths);
    tracked_hists_file->GetObject("hInnerAvReflectionResTimeVsCosTheta",hInnerAvReflectPaths);
    tracked_hists_file->GetObject("hMultipleEffectResTimeVsCosTheta",hMultipleEffectPaths);
    tracked_hists_file->GetObject("hAVPipesResTimeVsCosTheta",hAVPipesPaths);
    tracked_hists_file->GetObject("hAcrylicScatterResTimeVsCosTheta",hAcrylicPaths);
    tracked_hists_file->GetObject("OtherScatterResTimeVsCosTheta",hOtherScatterPaths);

    TH2F *hRegionCutReemissionResTimeVsCosTheta = (TH2F*)hReEmittedPaths->Clone();
    TH2F *hRegionCutPmtResTimeVsCosTheta = (TH2F*)hAllPaths->Clone();
    TH2F *hRegionCutNoisePaths = (TH2F*)hNoisePaths->Clone();
    TH2F *hRegionCutSingleScatterPaths = (TH2F*)hSingleScatterPaths->Clone();
    TH2F *hRegionCutOtherPaths = (TH2F*)hOtherPaths->Clone();
    TH2F *hRegionCutNoEffectPaths = (TH2F*)hNoEffectPaths->Clone();
    TH2F *hRegionCutNearReflectPaths = (TH2F*)hNearReflectPaths->Clone();
    TH2F *hRegionCutRopesPaths = (TH2F*)hRopesPaths->Clone();
    TH2F *hRegionCutPMTReflectionPaths = (TH2F*)hPMTReflectionPaths->Clone();
    TH2F *hRegionCutExtWaterScatterPaths = (TH2F*)hExtWaterScatterPaths->Clone();
    TH2F *hRegionCutInnerAvReflectPaths = (TH2F*)hInnerAvReflectPaths->Clone();
    TH2F *hRegionCutMultipleEffectPaths = (TH2F*)hMultipleEffectPaths->Clone();
    TH2F *hRegionCutAVPipesPaths = (TH2F*)hAVPipesPaths->Clone();
    TH2F *hRegionCutAcrylicPaths = (TH2F*)hAcrylicPaths->Clone();
    TH2F *hRegionCutOtherScatterPaths = (TH2F*)hOtherScatterPaths->Clone();

    hRegionCutReemissionResTimeVsCosTheta->SetName("hRegionSelectedReEmittedPaths");
    hRegionCutPmtResTimeVsCosTheta->SetName("hRegionSelectedAllPaths");
    hRegionCutNoisePaths->SetName("hRegionCutNoisePaths");
    hRegionCutSingleScatterPaths->SetName("hRegionCutSingleScatterPaths");
    hRegionCutOtherPaths->SetName("hRegionCutOtherPaths");
    hRegionCutNoEffectPaths->SetName("hRegionCutNoEffectPaths");
    hRegionCutNearReflectPaths->SetName("hRegionCutNearReflectPaths");
    hRegionCutRopesPaths->SetName("hRegionCutRopesPaths");
    hRegionCutPMTReflectionPaths->SetName("hRegionCutPMTReflectionPaths");
    hRegionCutExtWaterScatterPaths->SetName("hRegionCutExtWaterScatterPaths");
    hRegionCutInnerAvReflectPaths->SetName("hRegionCutInnerAvReflectPaths");
    hRegionCutMultipleEffectPaths->SetName("hRegionCutMultipleEffectPaths");
    hRegionCutAVPipesPaths->SetName("hRegionCutAVPipesPaths");
    hRegionCutAcrylicPaths->SetName("hRegionCutAcrylicPaths");
    hRegionCutOtherScatterPaths->SetName("hRegionCutOtherScatterPaths");

    int nBinsX = hReEmittedPaths->GetXaxis()->GetNbins();
    int nBinsY = hReEmittedPaths->GetYaxis()->GetNbins();

    for(int x=0; x<nBinsX+1; x++){
        double xBinCenter = hReEmittedPaths->GetXaxis()->GetBinCenter(x);
        for(int y=0; y<nBinsY+1; y++){
            double yBinCenter = hReEmittedPaths->GetYaxis()->GetBinCenter(y);
            double upperGradient = (y_b - y_a) / (x_b - x_a);
            double upperConstant = y_a - (upperGradient*x_a);
            double lowerGradient = (y_c - y_a) / (x_c - x_a);
            double lowerConstant = y_a - (lowerGradient*x_a);
            double rightGradient = 0;
            if(x_c != x_b) rightGradient = (y_c - y_b) / (x_c - x_b);
            double rightConstant = y_b - (rightGradient*x_b);
            if((yBinCenter <= upperGradient*xBinCenter + upperConstant) and (yBinCenter >= lowerGradient*xBinCenter + lowerConstant) and ((rightGradient < 0 and yBinCenter <= rightGradient*xBinCenter + rightConstant) or ((rightGradient > 0 and yBinCenter >= rightGradient*xBinCenter + rightConstant)) or (rightGradient == 0 and xBinCenter <= x_b))){
                hRegionCutReemissionResTimeVsCosTheta->SetBinContent(x,y,hReEmittedPaths->GetBinContent(x,y));
                hRegionCutPmtResTimeVsCosTheta->SetBinContent(x,y,hAllPaths->GetBinContent(x,y));
                hRegionCutNoisePaths->SetBinContent(x,y,hNoisePaths->GetBinContent(x,y));
                hRegionCutSingleScatterPaths->SetBinContent(x,y,hSingleScatterPaths->GetBinContent(x,y));
                hRegionCutOtherPaths->SetBinContent(x,y,hOtherPaths->GetBinContent(x,y));
                hRegionCutNoEffectPaths->SetBinContent(x,y,hNoEffectPaths->GetBinContent(x,y));
                hRegionCutNearReflectPaths->SetBinContent(x,y,hNearReflectPaths->GetBinContent(x,y));
                hRegionCutRopesPaths->SetBinContent(x,y,hRopesPaths->GetBinContent(x,y));
                hRegionCutPMTReflectionPaths->SetBinContent(x,y,hPMTReflectionPaths->GetBinContent(x,y));
                hRegionCutExtWaterScatterPaths->SetBinContent(x,y,hExtWaterScatterPaths->GetBinContent(x,y));
                hRegionCutInnerAvReflectPaths->SetBinContent(x,y,hInnerAvReflectPaths->GetBinContent(x,y));
                hRegionCutMultipleEffectPaths->SetBinContent(x,y,hMultipleEffectPaths->GetBinContent(x,y));
                hRegionCutAVPipesPaths->SetBinContent(x,y,hAVPipesPaths->GetBinContent(x,y));
                hRegionCutAcrylicPaths->SetBinContent(x,y,hAcrylicPaths->GetBinContent(x,y));
                hRegionCutOtherScatterPaths->SetBinContent(x,y,hOtherScatterPaths->GetBinContent(x,y));
            }
            else{
                hRegionCutReemissionResTimeVsCosTheta->SetBinContent(x,y,0);
                hRegionCutPmtResTimeVsCosTheta->SetBinContent(x,y,0);
                hRegionCutNoisePaths->SetBinContent(x,y,0);
                hRegionCutSingleScatterPaths->SetBinContent(x,y,0);
                hRegionCutOtherPaths->SetBinContent(x,y,0);
                hRegionCutNoEffectPaths->SetBinContent(x,y,0);
                hRegionCutNearReflectPaths->SetBinContent(x,y,0);
                hRegionCutRopesPaths->SetBinContent(x,y,0);
                hRegionCutPMTReflectionPaths->SetBinContent(x,y,0);
                hRegionCutExtWaterScatterPaths->SetBinContent(x,y,0);
                hRegionCutInnerAvReflectPaths->SetBinContent(x,y,0);
                hRegionCutMultipleEffectPaths->SetBinContent(x,y,0);
                hRegionCutAVPipesPaths->SetBinContent(x,y,0);
                hRegionCutAcrylicPaths->SetBinContent(x,y,0);
                hRegionCutOtherScatterPaths->SetBinContent(x,y,0);
            }
        }
    }

    //direct beam spot, max +/- 10 ns, -0.9 cos theta (straight through)

    TH2F *hDirectCutReemissionResTimeVsCosTheta = (TH2F*)hReEmittedPaths->Clone();
    TH2F *hDirectCutPmtResTimeVsCosTheta = (TH2F*)hAllPaths->Clone();
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

    hDirectCutReemissionResTimeVsCosTheta->SetName("hDirectCutReEmittedPaths");
    hDirectCutPmtResTimeVsCosTheta->SetName("hDirectCutAllPaths");
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

    double min_angle_direct_beam_spot = -0.9;

    for(int x=0; x<nBinsX+1; x++){
        double xBinCenter = hReEmittedPaths->GetXaxis()->GetBinCenter(x);
        if(xBinCenter <= min_angle_direct_beam_spot){
            for(int y=0; y<nBinsY+1; y++){
                double yBinCenter = hReEmittedPaths->GetYaxis()->GetBinCenter(y);
                if(yBinCenter >= min_time_direct_beam_spot and yBinCenter <= max_time_direct_beam_spot){
                    hDirectCutReemissionResTimeVsCosTheta->SetBinContent(x,y,hReEmittedPaths->GetBinContent(x,y));
                    hDirectCutPmtResTimeVsCosTheta->SetBinContent(x,y,hAllPaths->GetBinContent(x,y));
                    hDirectCutNoisePaths->SetBinContent(x,y,hNoisePaths->GetBinContent(x,y));
                    hDirectCutSingleScatterPaths->SetBinContent(x,y,hSingleScatterPaths->GetBinContent(x,y));
                    hDirectCutOtherPaths->SetBinContent(x,y,hOtherPaths->GetBinContent(x,y));
                    hDirectCutNoEffectPaths->SetBinContent(x,y,hNoEffectPaths->GetBinContent(x,y));
                    hDirectCutNearReflectPaths->SetBinContent(x,y,hNearReflectPaths->GetBinContent(x,y));
                    hDirectCutRopesPaths->SetBinContent(x,y,hRopesPaths->GetBinContent(x,y));
                    hDirectCutPMTReflectionPaths->SetBinContent(x,y,hPMTReflectionPaths->GetBinContent(x,y));
                    hDirectCutExtWaterScatterPaths->SetBinContent(x,y,hExtWaterScatterPaths->GetBinContent(x,y));
                    hDirectCutInnerAvReflectPaths->SetBinContent(x,y,hInnerAvReflectPaths->GetBinContent(x,y));
                    hDirectCutMultipleEffectPaths->SetBinContent(x,y,hMultipleEffectPaths->GetBinContent(x,y));
                    hDirectCutAVPipesPaths->SetBinContent(x,y,hAVPipesPaths->GetBinContent(x,y));
                    hDirectCutAcrylicPaths->SetBinContent(x,y,hAcrylicPaths->GetBinContent(x,y));
                    hDirectCutOtherScatterPaths->SetBinContent(x,y,hOtherScatterPaths->GetBinContent(x,y));
                }
                else{
                    hDirectCutReemissionResTimeVsCosTheta->SetBinContent(x,y,0);
                    hDirectCutPmtResTimeVsCosTheta->SetBinContent(x,y,0);
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
            }
        }
        else{
            for(int y=0; y<nBinsY+1; y++){
                hDirectCutReemissionResTimeVsCosTheta->SetBinContent(x,y,0);
                hDirectCutPmtResTimeVsCosTheta->SetBinContent(x,y,0);
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
        }
    }

    //reflected beam spot, max +/- 10 ns, 0.95 cos theta (straight through)

    TH2F *hReflectedCutReemissionResTimeVsCosTheta = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutPmtResTimeVsCosTheta = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutNoisePaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutSingleScatterPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutOtherPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutNoEffectPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutNearReflectPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutRopesPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutPMTReflectionPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutExtWaterScatterPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutInnerAvReflectPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutMultipleEffectPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutAVPipesPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutAcrylicPaths = (TH2F*)hAllPaths->Clone();
    TH2F *hReflectedCutOtherScatterPaths = (TH2F*)hAllPaths->Clone();

    hReflectedCutReemissionResTimeVsCosTheta->SetName("hReflectedCutReEmittedPaths");
    hReflectedCutPmtResTimeVsCosTheta->SetName("hReflectedCutAllPaths");
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

    double min_angle_reflected_beam_spot = 0.95;

    for(int x=0; x<nBinsX+1; x++){
        double xBinCenter = hReEmittedPaths->GetXaxis()->GetBinCenter(x);
        if(xBinCenter >= min_angle_reflected_beam_spot){
            for(int y=0; y<nBinsY+1; y++){
                double yBinCenter = hReEmittedPaths->GetYaxis()->GetBinCenter(y);
                if(yBinCenter >= min_time_reflected_beam_spot and yBinCenter <= max_time_reflected_beam_spot){
                    hReflectedCutReemissionResTimeVsCosTheta->SetBinContent(x,y,hReEmittedPaths->GetBinContent(x,y));
                    hReflectedCutPmtResTimeVsCosTheta->SetBinContent(x,y,hAllPaths->GetBinContent(x,y));
                    hReflectedCutNoisePaths->SetBinContent(x,y,hNoisePaths->GetBinContent(x,y));
                    hReflectedCutSingleScatterPaths->SetBinContent(x,y,hSingleScatterPaths->GetBinContent(x,y));
                    hReflectedCutOtherPaths->SetBinContent(x,y,hOtherPaths->GetBinContent(x,y));
                    hReflectedCutNoEffectPaths->SetBinContent(x,y,hNoEffectPaths->GetBinContent(x,y));
                    hReflectedCutNearReflectPaths->SetBinContent(x,y,hNearReflectPaths->GetBinContent(x,y));
                    hReflectedCutRopesPaths->SetBinContent(x,y,hRopesPaths->GetBinContent(x,y));
                    hReflectedCutPMTReflectionPaths->SetBinContent(x,y,hPMTReflectionPaths->GetBinContent(x,y));
                    hReflectedCutExtWaterScatterPaths->SetBinContent(x,y,hExtWaterScatterPaths->GetBinContent(x,y));
                    hReflectedCutInnerAvReflectPaths->SetBinContent(x,y,hInnerAvReflectPaths->GetBinContent(x,y));
                    hReflectedCutMultipleEffectPaths->SetBinContent(x,y,hMultipleEffectPaths->GetBinContent(x,y));
                    hReflectedCutAVPipesPaths->SetBinContent(x,y,hAVPipesPaths->GetBinContent(x,y));
                    hReflectedCutAcrylicPaths->SetBinContent(x,y,hAcrylicPaths->GetBinContent(x,y));
                    hReflectedCutOtherScatterPaths->SetBinContent(x,y,hOtherScatterPaths->GetBinContent(x,y));
                }
                else{
                    hReflectedCutReemissionResTimeVsCosTheta->SetBinContent(x,y,0);
                    hReflectedCutPmtResTimeVsCosTheta->SetBinContent(x,y,0);
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
        else{
            for(int y=0; y<nBinsY+1; y++){
                hReflectedCutReemissionResTimeVsCosTheta->SetBinContent(x,y,0);
                hReflectedCutPmtResTimeVsCosTheta->SetBinContent(x,y,0);
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

    //write to file

    std::string saveroot = "region_selected_hists_x_a" + std::to_string(x_a)  + "_x_b_" + std::to_string(x_b) + "_x_c_" + std::to_string(x_c) + "_y_a_" + std::to_string(y_a) + "_y_b_" + std::to_string(y_b) + "_y_c_" + std::to_string(y_c) + "_" + tracked_file;
    TFile *rootfile = new TFile(saveroot.c_str(),"RECREATE");
    rootfile->cd();
    hRegionCutReemissionResTimeVsCosTheta->Write();
    hRegionCutPmtResTimeVsCosTheta->Write();
    hRegionCutNoisePaths->Write();
    hRegionCutSingleScatterPaths->Write();
    hRegionCutOtherPaths->Write();
    hRegionCutNoEffectPaths->Write();
    hRegionCutNearReflectPaths->Write();
    hRegionCutRopesPaths->Write();
    hRegionCutPMTReflectionPaths->Write();
    hRegionCutExtWaterScatterPaths->Write();
    hRegionCutInnerAvReflectPaths->Write();
    hRegionCutMultipleEffectPaths->Write();
    hRegionCutAVPipesPaths->Write();
    hRegionCutAcrylicPaths->Write();
    hRegionCutOtherScatterPaths->Write();
    hDirectCutReemissionResTimeVsCosTheta->Write();
    hDirectCutPmtResTimeVsCosTheta->Write();
    hDirectCutNoisePaths->Write();
    hDirectCutSingleScatterPaths->Write();
    hDirectCutOtherPaths->Write();
    hDirectCutNoEffectPaths->Write();
    hDirectCutNearReflectPaths->Write();
    hDirectCutRopesPaths->Write();
    hDirectCutPMTReflectionPaths->Write();
    hDirectCutExtWaterScatterPaths->Write();
    hDirectCutInnerAvReflectPaths->Write();
    hDirectCutMultipleEffectPaths->Write();
    hDirectCutAVPipesPaths->Write();
    hDirectCutAcrylicPaths->Write();
    hDirectCutOtherScatterPaths->Write();
    hReflectedCutReemissionResTimeVsCosTheta->Write();
    hReflectedCutPmtResTimeVsCosTheta->Write();
    hReflectedCutNoisePaths->Write();
    hReflectedCutSingleScatterPaths->Write();
    hReflectedCutOtherPaths->Write();
    hReflectedCutNoEffectPaths->Write();
    hReflectedCutNearReflectPaths->Write();
    hReflectedCutRopesPaths->Write();
    hReflectedCutPMTReflectionPaths->Write();
    hReflectedCutExtWaterScatterPaths->Write();
    hReflectedCutInnerAvReflectPaths->Write();
    hReflectedCutMultipleEffectPaths->Write();
    hReflectedCutAVPipesPaths->Write();
    hReflectedCutAcrylicPaths->Write();
    hReflectedCutOtherScatterPaths->Write();
    rootfile->Write();
    rootfile->Close();

    return 0;
}