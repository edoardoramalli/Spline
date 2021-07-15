#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <thread>
#include <algorithm>
#include <iomanip>
#include <random>
//#include <dirent.h> // Linux



using namespace std;
using namespace std::chrono;



#include "Utilities.h"
#include "GlobalVariables.h"
#include "InputOutputSettings.h"
#include "BasisFunction.h"
#include "Spline.h"
#include "Indexes.h"
#include "dirent.h" // Windows
#include "CurveMatching.h"



int main() {

    // Reads the program settings from Settings.txt
    readProgramSettings();

    cout << "Curve Matching\n";

    // Displays the non-default Curve Matching settings and the graph settings
    displaySettings();

    // Runs the calculations and saves the results

    cout << "\nRunning";

    auto t1 = high_resolution_clock::now();

    CurveMatching curveMatching;

    curveMatching.run();

    auto t2 = high_resolution_clock::now();

    auto time = duration_cast<milliseconds>(t2-t1).count();
    for (int a=0; a<7; ++a) cout << "\b";
    cout << "Total calculation time: " << time << " ms\n\n";

    system("pause");

    return 0;

}