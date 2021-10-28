AMELLIE C++ scripts
===================

Contributors: Charlie

#### Requirements

- Compile with C++ 11
- Example compiler command: ```g++ -g -std=c++1y -o output.exe input.cpp `root-config --cflags --libs` -I$RATROOT/include/libpq -I$RATROOT/include -L$RATROOT/lib -lRATEvent_Linux```

## Usage

#### getTrackingInfo.cpp

This script generates a ROOT file containing histograms for the photon tracking information.

Requires three input arguments:
1. An AMELLIE simulation ROOT file, where the simulation saved all tracking information
2. The fibre which was simulated
3. The string "MC" or "raw" to describe the type of root file being used.

#### getRegions.cpp

This script takes the ROOT file generated in ```getTrackingInfo.cpp``` and selects the best triangular region for a given signal. 

Requires 6 input arguments:
1. The ROOT file containing the tracked histograms.
2. An integer specifying the number of bins to use in the optimisation. A larger number increases both precision and run time. The number of bins in the input histograms is the maximum by default, if you input a integer larger than that.
3. An integer representing a verbose flag (1 for true, 0 for false).
4. An integer representing a debug flag (1 for true, 0 for false).
5. An integer representing an extra info flag (1 for true, 0 for false). This will generate a few more histograms describing the evolution of the optimisation.
6. The required signal: use either ```reemitted```, ```scattered``` or ```attenuated```

The script will output a ROOT file containing histograms for these regions.

#### generateStats.cpp

This script has two functions which are accessed depending on the first input argument. The script can either generate stats for given regions or it can apply a region selection. The stats are produced in a txt file.

If wanting to generate stats, the script requires 4 input arguments:
1. The string ```generate_stats```
2. The ROOT file containing the region histograms (such as produced by ```getRegions.cpp``` or the other function in this script.) 
3. The ROOT file containing the tracked histograms, as produced by ```getTrackingInfo.cpp```.
4. The required signal: use either ```reemitted```, ```scattered``` or ```attenuated```
5. The string "MC" or "raw" to describe the type of root file being used.

If wanting to apply a region selection, the script requires 12 input arguments:
1. The string ```apply_region```
2. The ROOT file containing the tracked histograms, as produced by ```getTrackingInfo.cpp```
3. A number corresponding to x_a - The cos theta value of the left most triangular point
4. A number corresponding to x_b - The cos theta value of the upper right most triangular point
5. A number corresponding to x_c - The cos theta value of the lower right most triangular point
6. A number corresponding to y_a - The residual time value of the left most triangular point
7. A number corresponding to y_b - The residual time value of the upper right most triangular point
8. A number corresponding to y_c - The residual time value of the lower right most triangular point
9. The minimum residual time for the direct beam spot box cut
10. The maximum residual time for the direct beam spot box cut
11. The minimum residual time for the reflected beam spot box cut
12. The maximum residual time for the reflected beam spot box cut
13. The string "MC" or "raw" to describe the type of root file being used.
