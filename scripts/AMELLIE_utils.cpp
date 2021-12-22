#include "AMELLIE_utils.hpp"


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/**
 * @brief Construct a new triangle::triangle object
 * 
 * @param x_a 
 * @param x_b 
 * @param x_c 
 * @param y_a 
 * @param y_b 
 * @param y_c 
 */
triangle::triangle(double x_a, double x_b, double x_c, double y_a, double y_b, double y_c) {
    points[0] = x_a; points[1] = x_b; points[2] = x_c;
    points[3] = y_a; points[4] = y_b; points[5] = y_c;
    // Twice the area of the triangle (twice to save on computation later)
    Area2 = abs(points[0] * (points[4] - points[5]) + points[1] * (points[5] - points[3])
                + points[2] * (points[3] - points[4]));
}

// member functions

// Get/Set triangle vertex coordinates
double& triangle::X_a() {return points[0];}
double& triangle::X_b() {return points[1];}
double& triangle::X_c() {return points[2];}
double& triangle::Y_a() {return points[3];}
double& triangle::Y_b() {return points[4];}
double& triangle::Y_c() {return points[5];}

double triangle::Area() {return 0.5 * Area2;}

/**
 * @brief Checks if point is inside triangle by using the sum of triangles method:
 * if the sum of the area of the three triangles formed by the point and two vertices of the triangle add up to
 * the area of the original triangle, the point is inside. It is outside if the sum is larger. 
 * 
 * @param point_x 
 * @param point_y 
 * @param lim_count true (default) = point on triangle edge counts as inside, false = point on edge counts as outside.
 * @return true 
 * @return false 
 */
bool triangle::check_point_inside_triangle(const double point_x, const double point_y, const bool lim_count) {
    if (Area2 == 0.0) {
        // std::cout << "Triangle area is zero." << std::endl;
        if (point_x == points[0] and point_y == points[3]) {
            return lim_count;  // still might be on the point of the triangle
        } else {
            return false;
        }
    }

    // Calculate (twice) the area of each triangle formed by the point, and check if it is zero successively
    double A_a = abs(point_x * (points[4] - points[5]) + points[1] * (points[5] - point_y)
                    + points[2] * (point_y - points[4]));
    if (A_a == 0.0) {return lim_count;} // if area is zero, point is on edge of triangle.

    double A_b = abs(points[0] * (point_y - points[5]) + point_x * (points[5] - points[3])
                    + points[2] * (points[3] - point_y));
    if (A_b == 0.0) {return lim_count;} // if area is zero, point is on edge of triangle.

    double A_c = abs(points[0] * (points[4] - point_y) + points[1] * (point_y - points[3])
                    + point_x * (points[3] - points[4]));
    if (A_c == 0.0) {return lim_count;} // if area is zero, point is on edge of triangle.

    // Check sum of areas of point triangles adds up to original triangle aread (tolerance of 0.1%).
    // If they equal, the point is inside triangle. If the sum is larger, the point is outside.
    double frac_diff = ((A_a + A_b + A_c) / Area2) - 1.0;
    if (frac_diff <= 0.001) {
        return true;
    } else {
        return false;
    }
}


/* ~~~~~~~~~ operator overload ~~~~~~~~~ */
double& triangle::operator [] (int i) {
    assert((i < 6 && i >= 0) && "Index out of range");
    return points[i];
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HISTLIST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/**
 * @brief Construct a new Hist List:: Hist List object.
 * Default constructor
 * 
 */
HistList::HistList() {length = 0;}

/**
 * @brief Construct a new Hist List:: Hist List object.
 * Creates vector of 2D Tracking Hists in root file, and clones them to
 * 3 new vectors of 2D hists with modified names.
 * 
 * @param tracking_file Root file containing hists.
 * @param name_list List of histogram names the user wishes to read in. Default list defined
 * in header file. (FIXME: can add option to just read in all 2D hists)
 */
HistList::HistList(std::string tracking_file, std::vector<std::string> name_list) {
    Read_File(tracking_file, name_list);
}

/**
 * @brief Basically constructs a new HistList object.
 * Creates vector of 2D Tracking Hists in root file, and clones them to
 * 3 new vectors of 2D hists with modified names.
 * 
 * @param tracking_file Root file containing hists.
 * @param name_list List of histogram names the user wishes to read in. Default list defined
 * in header file. (FIXME: can add option to just read in all 2D hists)
 */
void HistList::Read_File(std::string tracking_file, std::vector<std::string> name_list) {
    // Read in root file
    TFile* fin = new TFile(tracking_file.c_str());
    if (!fin->IsOpen()) {
        std::cout << "Cannot open input file " << tracking_file << std::endl;
        exit(1);
    }
    length = name_list.size();
    
    // Iterate through list of objects in root file
    TList* list = fin->GetListOfKeys() ;
    if (!list) {std::cout << "No keys found in file\n" << std::endl; exit(1);}
    TIter next(list);
    TObject* obj;
    TKey* key;
    int j = 0;
    std::string name;
    int order[name_list.size()];
    for (unsigned int i = 0; i < length; ++i) {
        order[i] = -1;
    }
    std::vector<TH2F*> temp_HistList;

    // Go through list of histograms and add them to temp list if they are included in name list
    while((key = (TKey*)next())){
        obj = key->ReadObj() ;
        if(obj->InheritsFrom(TH2::Class())){
            // Check which histogram in file matches name
            name = (std::string)(((TH2F*)obj)->GetName());
            for (unsigned int i = 0; i < length; ++i) {
                if (name == name_list.at(i)) {
                    // Add tracking histograms to temporary list
                    temp_HistList.push_back((TH2F*)obj);
                    order[i] = j;
                    ++j;
                    break;
                }
            }
        }
    }

    // Reorder histograms, and clone new ones in that order
    for (unsigned int i = 0; i < length; ++i) {
        if (order[i] == -1) {
            std::cout << "Histogram '" << name_list.at(i) << "' not found." << std::endl;
            exit(1);
        } else {
            // Add tracking histograms to list
            tracking_hists.push_back(temp_HistList.at(order[i]));

            // Clone tracking hists to other lists
            name = name_list.at(i).erase(0,1); // removes first character ("h" in this case)
            region_hists.push_back((TH2F*)(tracking_hists.at(i))->Clone());
            region_hists.at(i)->SetName(((std::string)"hRegion" + name).c_str());
            direct_hists.push_back((TH2F*)(tracking_hists.at(i))->Clone());
            direct_hists.at(i)->SetName(((std::string)"hDirect" + name).c_str());
            reflected_hists.push_back((TH2F*)(tracking_hists.at(i))->Clone());
            reflected_hists.at(i)->SetName(((std::string)"hReflected" + name).c_str());
        }
    }
}

// return suitable histograms / histogram names
unsigned int HistList::len() {return length;}
std::vector<TH2F*>& HistList::Tracking_Hists() {return tracking_hists;}
std::vector<TH2F*>& HistList::Region_Hists() {return region_hists;}
std::vector<TH2F*>& HistList::Direct_Hists() {return direct_hists;}
std::vector<TH2F*>& HistList::Reflected_Hists() {return reflected_hists;}

/**
 * @brief Write all  new (Region, Direct, Reflected) hists to file
 * 
 */
void HistList::Write() {
    for (unsigned int i = 0; i < length; ++i) {
        region_hists.at(i)->Write();
        direct_hists.at(i)->Write();
        reflected_hists.at(i)->Write();
    }
}