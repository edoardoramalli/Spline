class Indexes {

public:

    /* Name of the _exp.txt and _mod.txt files with the data for the splines */
    string fileName;

    /* Number of models to be compared with the experimental data */
    int numberOfModels;

    /* Names of the models */
    vector<string> modelNames;

    /* Scores, before the shift. scores[i] contains the data relative to the
    comparison of modelNames[i] and the experimental data. scores[i] is equal to
    -1 in case it is not possible to calculate score i, or in case the
    calculated score is nan */
    vector<double> scores;

    /* Scores, after the shift. scores_shift[i] contains the data relative to
    the comparison of modelNames[i] and the experimental data. scores_shift[i]
    is equal to -1 in case it is not possible to calculate score i, or in case
    the calculated score is nan */
    vector<double> scores_shift;

    /* Standard deviations of the scores, before the shift. stdDev[i] contains
    the data relative to the comparison of modelNames[i] and the experimental
    data. stdDev[i] is equal to -1 in case it is not possible to calculate the
    corresponding score, or it is nan */
    vector<double> stdDev;

    /* Standard deviations of the scores, after the shift. stdDev_shift[i]
    contains the data relative to the comparison of modelNames[i] and the
    experimental data. stdDev_shift[i] is equal to -1 in case it is not possible
    to calculate the corresponding score, or it is nan */
    vector<double> stdDev_shift;

    /* Individual indexes for each variation of the experimental data obtained
    with the bootstrapping procedure. allIndexes[i][j][k]: i: variations of the
    experimental data; j: models; k: d0L2, d1L2, d0Pe, d1Pe, shift (or shifts if
    all different) */
    vector<vector<vector<double>>> allIndexes;

    ////////////////////////////////////////////////////////////////////////////

    /* Executes all the operations for the comparison of the experimental data
    with the models */
    void solve(const string& folderPath,
               const string& folderName,
               const string& fileName,
               bool calculateShift);

    ////////////////////////////////////////////////////////////////////////////

private:

    /* Equivalent of 'scores' for every bootstrap variation */
    vector<vector<double>> bootstrapScores;

    /* Equivalent of scores_shift for every bootstrap variation */
    vector<vector<double>> bootstrapScores_shift;

    /* x-coordinates and y-coordinates of the experimental data and the models.
    inputData[0] and inputData[1] contain, respectively, the x-coordinates and
    the y-coordinates of the experimental data, inputData[i] and inputData[i+1]
    contain, respectively, the x-coordinates and the y-coordinates of
    modelNames[i/2-1] */
    vector<vector<double>> inputData;

    /* Splines for the experimental data */
    vector<Spline> splinesExp;

    /* Splines for the models. splinesMod[i] contains the spline for
    modelNames[i] */
    vector<Spline> splinesMod;

    /* Splines for the models, before being extended */
    vector<Spline> splinesMod_original;

    /* Path to the folder containing the input data */
    string folderPath;

    /* Name of the folder containing the input data */
    string folderName;

    /* Relative errors for each experimental data point */
    vector<double> relativeErrors;

    /* Specifies whether to calculate the indexes for the shifted splines */
    bool calculateShift;

    /* Amounts to be added to the abscissae of the models to find the shifted
    abscissae. shiftAmounts[i][j] refers to splinesExp[i] and modelNames[j] */
    vector<vector<vector<double>>> shiftAmounts;

    /* Positions in shiftAmounts, indicating the shifts that maximise the mean
    of d0L2, d1L2, d0Pe and d1Pe for the models. maxIndexes[i][j] refers to
    shiftAmounts[i][j] */
    vector<vector<int>> maxIndexes;

    /* Index in the splinesMod vector of the model being calculated */
    int indexOfModelBeingCalculated;

    /* Height of the trapezoids used used during integration, divided by 2 */
    double halfHeight;

    /* Powers of the abscissae of the sides of the trapezoids used for the
    numerical calculation of the indexes, distributed evenly between the
    extremes of the experimental data */
    vector<vector<double>> xPowers;

    /* Ordinates of the experimental data spline, with yD0[i] corresponding to
    xPowers[i][1] */
    vector<double> yD0;

    /* Ordinates of the first derivative of the experimental data spline, with
    yD1[i] corresponding to xPowers[i][1] */
    vector<double> yD1;

    /* Scores without shifts. Sorted from lowest to highest. Is equal to -1 in
    case it is not possible to calculate a score */
    vector<double> scoresSorted;

    /* Standard deviations without shifts. Sorted from lowest to highest. Is
    equal to -1 in case it is not possible to calculate a score */
    vector<double> stdDevSorted;

    /* Model names, sorted as the elements of scoresSorted */
    vector<string> modelNamesSorted;

    /* Line types in the .R graphs, sorted as the elements of scoresSorted */
    vector<int> scoresSortedLineTypes;

    /* Color indexes in the .R graphs, sorted as the elements of scoresSorted */
    vector<int> scoresSortedColorIndexes;

    /* Scores with shifts. Sorted from lowest to highest. Is equal to -1 in case
    it is not possible to calculate a score */
    vector<double> scoresSorted_shift;

    /* Standard deviations with shifts. Sorted from lowest to highest. Is equal
    to -1 in case it is not possible to calculate a score */
    vector<double> stdDevSorted_shift;

    /* Model names, sorted as the elements of scoresSorted_shift */
    vector<string> modelNamesSorted_shift;

    /* Line types in the .R graphs, sorted as the elements of scoresSorted_shift
    */
    vector<int> scoresSortedLineType_shift;

    /* Color indexes in the .R graphs, sorted as the elements of
    scoresSorted_shift */
    vector<int> scoresSortedColorIndexes_shift;

    /* Index in the splinesExp vector of the spline being considered */
    int splineExpIndex;

    /* Index in the splinesExp vector of the best fit of the experimental data
    */
    int indexBestSplineExp;

    /* Index in the splinesExp vector of the best fit of the experimental data
    */
    int newindexBestSplineExp;

    /* Values of the d0L2, d1L2, d0Pe and d1Pe indexes, in this order, obtained
    before any shift. indexes[i][j] refers to splinesExp[i] and contains the
    data relative to the comparison of modelNames[j] and the experimental data.
    Is equal to -1 in case it is not possible to calculate an index */
    vector<vector<vector<double>>> indexes;

    /* Values of the d0L2_shift, d1L2_shift, d0Pe_shift, d1Pe_shift indexes, in
    this order, corresponding to their optimal values found while shifting.
    indexes_shift[i][j] refers to splinesExp[i] and contains the data relative
    to the comparison of modelNames[j] and the experimental data. Is equal to -1
    in case it is not possible to calculate an index */
    vector<vector<vector<double>>> indexes_shift;

    /* Values of the shift indexes, corresponding to the normalized shift
    amounts which maximize the means of d0L2, d1L2, d0Pe and d1Pe for each
    model. shifts[i][j] refers to splinesExp[i] and contains the data relative
    to the comparison of modelNames[j] and the experimental data. Is equal to -1
    in case it is not possible to calculate a shift */
    vector<vector<double>> shifts;

    /* Values of the shift indexes, corresponding to the normalized shift
    amounts which maximise each of d0L2, d1L2, d0Pe and d1Pe for each model.
    shifts[i][j] refers to splinesExp[i] and contains the data relative to the
    comparison of modelNames[j] and the experimental data. Is equal to -1 in
    case it is not possible to calculate a shift */
    vector<vector<vector<double>>> shifts_allDifferent;

    /* Norm of the experimental data, calculated between the extremes of the
    experimental data */
    double normD0Exp;

    /* Norm of the first derivative of the experimental data, calculated between
    the extremes of the experimental data */
    double normD1Exp;

    /* Norm of a model, calculated between the extremes of the experimental data
    */
    double normD0Mod;

    /* Norm of the first derivative of a model, calculated between the extremes
    of the experimental data */
    double normD1Mod;

    /* Specifies whether for a model it is possible to calculate the indexes for
    the shifted model */
    vector<vector<bool>> possibleToCalculateShifts;

    /* Integral of the absolute value of the first derivative of the
    experimental data */
    double integralOfAbsoluteValueD1Exp;

    /* Index of the bootstrap variation currently being considered */
    int bootstrapIndex;

    /* Used by calculateShiftAmounts */
    double maxExtrapolatedLength;

    /* Used by calculateShiftAmounts */
    double shiftAroundMaximum;

    /* Used by calculateShiftAmounts */
    double distanceBetweenShiftedPoints;

    ////////////////////////////////////////////////////////////////////////////

    /* Reads the data relative to the experiment (x, y, relative errors) and the
    models (names, x, y). Sorts the data in ascending order relative to the
    abscissae. If two or more abscissae are the same, replaces them with a
    single abscissa with the mean of the corresponding ordinates as ordinate,
    and with the mean of the corresponding relative errors as relative error.
    The data is read from folderPath/fileName_exp.txt and from
    folderPath/fileName_mod.txt */
    void readData();

    /* Calculates the splines for the experimental data */
    void calculateSplines();

    /* Calculates the indexes for each model */
    void calculateIndexes();

    /* Calculates the scores from the indexes */
    void calculateScores();

    /* Saves the data necessary to plot the splines to .txt files */
    void saveGraphData();
    
    /*Calculate the summed squared error to select the best spline for the 
    experimental data*/

    double Stat_SSE(vector <double> b, vector<double> c);

};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



void Indexes::solve(const string& FolderPath,const string& FolderName,const string& FileName,bool CalculateShift) {

    folderPath = FolderPath;
    folderName = FolderName;
    fileName = FileName;
    calculateShift = CalculateShift;

    // Processes the data in the input files
    this->readData();


    // Calculates the splines for the experimental data
    this->calculateSplines();

    // For each knot arrangement in the experimental data, calculates the
    // indexes for each model
    for (int b=0; b<splinesExp.size(); ++b) {
        splineExpIndex = b;
        this->calculateIndexes();
    }

    // Calculates the scores from the indexes
    this->calculateScores();

    

    // Saves the data necessary to plot the splines to .txt files.
    // Used to interface with Alice's work
    this->saveGraphData();

}



void Indexes::readData() {

    // Copies the data from folderPath/name_exp.txt to dataStringExp

    vector<vector<string>> dataStringExp;

    string fileNameWithPathExp = folderPath + fileName + "_exp.txt";
    const char* fileNameWithPathCharExp = fileNameWithPathExp.c_str();
    string lineExp;
    string elementExp;
    vector<string> emptyVectorExp;
    int rowIndexExp = 0;
    ifstream myfileExp(fileNameWithPathCharExp);
    while (getline(myfileExp,lineExp)) {
        dataStringExp.push_back(emptyVectorExp);
        stringstream lineStream(lineExp);
        while (lineStream >> elementExp)
            dataStringExp[rowIndexExp].push_back(elementExp);
        ++rowIndexExp;
    }
    myfileExp.close();

    // Converts the input data from string to double and saves it to
    // dataDoubleExp, with the data saved in reverse order relative to the order
    // in the input file (to speed up processing)

    auto dataDoubleExp = vector<vector<double>>(3);

    for (int a=dataStringExp.size()-1; a>0; --a)
        if (dataStringExp[a][0] != "NA") {
            dataDoubleExp[0].push_back(stod(dataStringExp[a][0]));
            dataDoubleExp[1].push_back(stod(dataStringExp[a][1]));
            if (dataStringExp[0].size() == 3)
                dataDoubleExp[2].push_back(stod(dataStringExp[a][2]));
            else
                dataDoubleExp[2].push_back(defaultRelativeError);
        }

    // Fills dataSortedExp, corresponding to dataDoubleExp with the data sorted
    // in ascending order with respect to the abscissae

    auto dataSortedExp = vector<vector<double>>(3);

    while (dataDoubleExp[0].size() > 0)
        for (int a=dataDoubleExp[0].size()-1; a>-1; --a) {
            double value = dataDoubleExp[0][a];
            int index = a;
            for (int b=dataDoubleExp[0].size()-1; b>-1; --b)
                if (dataDoubleExp[0][b] < value) {
                    value = dataDoubleExp[0][b];
                    index = b;
                }
            dataSortedExp[0].push_back(dataDoubleExp[0][index]);
            dataSortedExp[1].push_back(dataDoubleExp[1][index]);
            dataSortedExp[2].push_back(dataDoubleExp[2][index]);
            dataDoubleExp[0].erase(dataDoubleExp[0].begin()+index);
            dataDoubleExp[1].erase(dataDoubleExp[1].begin()+index);
            dataDoubleExp[2].erase(dataDoubleExp[2].begin()+index);
        }

    

    // If two or more abscissae are the same, replaces them with a single
    // abscissa with the mean of the corresponding ordinates as ordinate, and
    // (for experimental data) with the mean of the corresponding relative
    // errors as relative error. Copies the resulting abscissae and ordinates to
    // inputData and the resulting relative errors to allRelativeErrors

    auto allInputData = vector<vector<double>>(dataSortedExp.size());
    vector<double> allRelativeErrors;

    allInputData[0].push_back(dataSortedExp[0][0]);
    allInputData[1].push_back(dataSortedExp[1][0]);
    allRelativeErrors.push_back(dataSortedExp[2][0]);

    int nExp = 1;

    for (int a=1; a<dataSortedExp[0].size(); ++a) {
        if (dataSortedExp[0][a] != allInputData[0].back()) {
            allInputData[0].push_back(dataSortedExp[0][a]);
            allInputData[1].push_back(dataSortedExp[1][a]);
            allRelativeErrors.push_back(dataSortedExp[2][a]);
            nExp = 1;
        }
        else {
            allInputData[1].back() =
                allInputData[1].back()*(double)nExp + dataSortedExp[1][a];
            allRelativeErrors.back() =
                allRelativeErrors.back()*(double)nExp + dataSortedExp[2][a];
            ++nExp;
            allInputData[1].back() /= (double)nExp;
            allRelativeErrors.back() /= (double)nExp;
        }
    }


    // Obtains inputData and relativeErrors by removing excess points with
    // ordinate 0 on the left and on the right of experimental data and models
    
    inputData = vector<vector<double>>(allInputData.size());

    for (int a=0; a<allInputData.size(); ++a)
        if (allInputData[a].size() > 1) {

            if (a == 0) { // qui sta selzionando le ascisse sperimentali

                double maxOrdinate = allInputData[a+1][0];
                double minOrdinate = allInputData[a+1][0];
                for (int b=1; b<allInputData[a].size(); ++b) {
                    if (allInputData[a+1][b] > maxOrdinate)
                        maxOrdinate = allInputData[a+1][b];
                    if (allInputData[a+1][b] < minOrdinate)
                        minOrdinate = allInputData[a+1][b];
                }
                double height = maxOrdinate - minOrdinate;
                double minHeight =
                    fractionOfOrdinateRangeForAsymptoteIdentification * height;
                int indexOfFirstPoint = 0;
                for (int b=1; b<allInputData[a].size(); ++b)
                    if (fabs(allInputData[a+1][b]-allInputData[a+1][b-1]) >
                        minHeight) {
                        indexOfFirstPoint = b-1;
                        break;
                    }
                int indexOfLastPoint = allInputData[a].size()-1;
                for (int b=allInputData[a].size()-2; b>-1; --b)
                    if (fabs(allInputData[a+1][b]-allInputData[a+1][b+1]) >
                        minHeight) {
                        indexOfLastPoint = b+1;
                        break;
                    }

                inputData[a] =
                    vector<double>(indexOfLastPoint-indexOfFirstPoint+1);
                inputData[a+1] =
                    vector<double>(indexOfLastPoint-indexOfFirstPoint+1);

                for (int b=0; b<inputData[a].size(); ++b)
                    inputData[a][b] = allInputData[a][indexOfFirstPoint+b];

                for (int b=0; b<inputData[a].size(); ++b)
                    inputData[a+1][b] = allInputData[a+1][indexOfFirstPoint+b];

                for (int b=indexOfFirstPoint; b<=indexOfLastPoint; ++b)
                    relativeErrors.push_back(allRelativeErrors[b]);

            }

            if (a > 0 && a%2 == 0) {

                inputData[a] = allInputData[a];
                inputData[a+1] = allInputData[a+1];

            }

        }

}



void Indexes::calculateSplines() {

    if (inputData[0].size() < 3)
        splinesExp = vector<Spline>(1);
    else if (inputData[0].size() < 5)
        splinesExp = vector<Spline>(2);
    else
        splinesExp = vector<Spline>(3);

    splinesExp[0].solve(inputData[0],inputData[1],0,0);
    splinesExp[0].removeAsymptotes();

    if (splinesExp.size() > 1) {
        splinesExp[1].solve(inputData[0],inputData[1],0,2);
        splinesExp[1].removeAsymptotes();
    }

    if (splinesExp.size() > 2) {
        splinesExp[2].solve(inputData[0],inputData[1],0,5);
        splinesExp[2].removeAsymptotes();
    }



}



void Indexes::calculateIndexes() {

    int i = splineExpIndex;

    double minLengthInCommon = splinesExp[i].xRange*fractionOfExpRangeForModelExtrapolation ;

    
    // DA VEDERE
    //maxExtrapolatedLength = splinesExp[i].xRange - minLengthInCommon

    // Normalizes the spline coefficients of the experimental data and of the
    // first derivative of the experimental data
    splinesExp[i].normalizeCoefficients(
                                    -splinesExp[i].yD0Min,
                                    splinesExp[i].yD0Max-splinesExp[i].yD0Min,
                                    splinesExp[i].yD1MaxAbs);

    

}



void Indexes::calculateScores() {

    // Selects the best spline of splinesExp

    indexBestSplineExp = 0;

    vector<double> ySpl_0;
    vector<double> ySpl_1;
    vector<double> ySpl_2;
    
    for (int i=0; i<inputData[0].size();i++)
        ySpl_0.push_back(splinesExp[0].D0(inputData[0][i]));
    
    for (int i=0; i<inputData[0].size();i++)
        ySpl_1.push_back(splinesExp[1].D0(inputData[0][i]));

    for (int i=0; i<inputData[0].size();i++)
        ySpl_2.push_back(splinesExp[2].D0(inputData[0][i]));
    
    double SSE_0 = Stat_SSE(inputData[1],ySpl_0);

    double SSE_1 = Stat_SSE(inputData[1],ySpl_1);

    double SSE_2 = Stat_SSE(inputData[1],ySpl_2);

    if (SSE_0<=SSE_1 && SSE_0<=SSE_2){
        indexBestSplineExp = 0;
    }
    if (SSE_1<=SSE_0 && SSE_1<=SSE_2){
        indexBestSplineExp = 1;
    }
    if (SSE_2<=SSE_0 && SSE_2<=SSE_1){
        indexBestSplineExp = 2;
    }
      
    if(inputData[0].back() != splinesExp[indexBestSplineExp].knots.back() || inputData[0][0] != splinesExp[indexBestSplineExp].knots[0]){

        if(indexBestSplineExp == 0){
            if (SSE_1<=SSE_2){
                newindexBestSplineExp = 1;
            }
            
            if (SSE_2<=SSE_0){
                newindexBestSplineExp = 2;
            }
        }
        if(indexBestSplineExp == 1){
            if (SSE_0<=SSE_2){
                newindexBestSplineExp = 0;
            }
            
            if (SSE_2<=SSE_0 ){
                newindexBestSplineExp = 2;
            }
        }
        if(indexBestSplineExp == 2){
            if (SSE_0<=SSE_1){
                newindexBestSplineExp = 0;
            }
            
            if (SSE_1<=SSE_0 ){
                newindexBestSplineExp = 1;
            }
        }
    }

    else 
        newindexBestSplineExp = indexBestSplineExp;


    indexBestSplineExp = newindexBestSplineExp;    
}


void Indexes::saveGraphData() {

    if (graphs == false) return;

    int j = indexBestSplineExp;
    int p;
    int q;

    auto x = vector<double>(graphPoints);
    auto x_1 = vector<double>(graphPoints);
    auto x_2 = vector<double>(graphPoints);
    vector<double> xShift;
    vector<vector<double>> xShiftPowers;
    if (graphsD0Shift == true || graphsD1Shift == true) {
        xShift = vector<double>(graphPoints);
        xShiftPowers = vector<vector<double>>(graphPoints,vector<double>(m,1));
    }

    vector<double> y_D0, y_D1, y_D0_Shift, y_D1_Shift, nodi;
    vector<double> y_D0_1, nodi_1;
    vector<double> y_D0_2, nodi_2;
    if (graphsD0 == true) y_D0 = vector<double>(graphPoints);
    if (graphsD1 == true) y_D1 = vector<double>(graphPoints);
    if (graphsD0Shift == true) y_D0_Shift = vector<double>(graphPoints);
    if (graphsD1Shift == true) y_D1_Shift = vector<double>(graphPoints);
    if (knotsToSave == true) nodi = vector<double>(graphPoints);
    if (allSplineExp == true) y_D0_1 = vector<double>(graphPoints), nodi_1 = vector<double>(graphPoints), 
                                y_D0_2 = vector<double>(graphPoints), nodi_2 = vector<double>(graphPoints);


    if (splinesExp[j].possibleToCalculateSpline == true) {

        // Selects the points on the x-axis

        double distance =
            (splinesExp[j].knots.back()-splinesExp[j].knots[0]) /
            (double)(graphPoints);

        for (int b=0; b<graphPoints; ++b)
            x[b] = splinesExp[j].knots[0]+(double)b*distance;
        x.back() = splinesExp[j].knots.back();

        // Calculates the ordinates

        if (graphsD0 == true)
            for (int b=0; b<graphPoints; ++b)
                y_D0[b] = splinesExp[j].D0(x[b]);

        if (graphsD1 == true)
            for (int b=0; b<graphPoints; ++b)
                y_D1[b] = splinesExp[j].D1(x[b]);

        if (graphsD0Shift == true) {
            if (graphsD0 == true)
                y_D0_Shift = y_D0;
            else
                for (int b=0; b<graphPoints; ++b)
                    y_D0_Shift[b] = splinesExp[j].D0(x[b]);
        }

        if (graphsD1Shift == true) {
            if (graphsD1 == true)
                y_D1_Shift = y_D1;
            else
                for (int b=0; b<graphPoints; ++b)
                    y_D1_Shift[b] = splinesExp[j].D1(x[b]);
        }
        if (knotsToSave == true)
            for (int b=0; b<splinesExp[j].knots.size(); ++b)
                nodi[b] = splinesExp[j].knots[b];

        if (allSplineExp == true){

                if (j == 0){
                    p = 1;
                    q = 2;
                }
                else if (j == 1){
                    p = 0;
                    q = 2;
                }
                else if  (j == 2){
                    p = 0;
                    q = 1;
                }
                

            double distance_1 = (splinesExp[p].knots.back()-splinesExp[p].knots[0]) / (double)(graphPoints-1);

            for (int b=0; b<graphPoints-1; ++b)
                x_1[b] = splinesExp[p].knots[0]+(double)b*distance_1;
            x_1.back() = splinesExp[p].knots.back();

            double distance_2 = (splinesExp[q].knots.back()-splinesExp[q].knots[0]) / (double)(graphPoints-1);

            for (int b=0; b<graphPoints-1; ++b)
                x_2[b] = splinesExp[q].knots[0]+(double)b*distance_2;
            x_2.back() = splinesExp[q].knots.back();

            for (int b=0; b<graphPoints; ++b)
                y_D0_1[b] = splinesExp[p].D0(x_1[b]);

            for (int b=0; b<graphPoints; ++b)
                y_D0_2[b] = splinesExp[q].D0(x_2[b]);
            
            for (int b=0; b<splinesExp[p].knots.size(); ++b)
                nodi_1[b] = splinesExp[p].knots[b];
            
            for (int b=0; b<splinesExp[q].knots.size(); ++b)
                nodi_2[b] = splinesExp[q].knots[b];

        } // End of allSpline 

        // Saves the results to .txt files

        string pathAndPartialName =
            "./Graph data/" + fileName + "_";

        if (graphsD0 == true) {

            string pathAndNameString = pathAndPartialName + "Exp_D0.txt";

            const char* pathAndName = pathAndNameString.c_str();
            ofstream script;
            script.open(pathAndName,std::ios::app);

            script << "x_Exp\t" << fileName << "_Exp\n";
            for (int b=0; b<graphPoints; ++b)
                script << x[b] << "\t" << y_D0[b] << "\n";

            script.close();

        }

        if (saveExpData == true) {

            string pathAndNameString = pathAndPartialName + "Exp_Data.txt";

            const char* pathAndName = pathAndNameString.c_str();
            ofstream script;
            script.open(pathAndName,std::ios::app);

            script << "x_Exp\t" << "y_Exp\n";
            for (int b=0; b<inputData[0].size(); ++b)
                script << inputData[0][b] << "\t" << inputData[1][b]  << "\n";

            script.close();

        }

        if (graphsD1 == true) {

            string pathAndNameString = pathAndPartialName + "Exp_D1.txt";

            const char* pathAndName = pathAndNameString.c_str();
            ofstream script;
            script.open(pathAndName,std::ios::app);

            script << "x_Exp\t" << fileName << "_Exp\n";
            for (int b=0; b<graphPoints; ++b)
                script << x[b] << "\t" << y_D1[b] << "\n";

            script.close();

        }

        if (graphsD0Shift == true) {

            string pathAndNameString = pathAndPartialName + "Exp_D0_shift.txt";

            const char* pathAndName = pathAndNameString.c_str();
            ofstream script;
            script.open(pathAndName,std::ios::app);

            script << "x_Exp\t" << fileName << "_Exp\n";
            for (int b=0; b<graphPoints; ++b)
                script << x[b] << "\t" << y_D0_Shift[b] << "\n";

            script.close();

        }

        if (graphsD1Shift == true) {

            string pathAndNameString = pathAndPartialName + "Exp_D1_shift.txt";

            const char* pathAndName = pathAndNameString.c_str();
            ofstream script;
            script.open(pathAndName,std::ios::app);

            script << "x_Exp\t" << fileName << "_Exp\n";
            for (int b=0; b<graphPoints; ++b)
                script << x[b] << "\t" << y_D1_Shift[b] << "\n";

            script.close();

        }

        if (knotsToSave == true) {

            string pathAndNameString = pathAndPartialName + "knots.txt";

            const char* pathAndName = pathAndNameString.c_str();
            ofstream script;
            script.open(pathAndName,std::ios::app);

            script << "knots\t"<<"indexBestSplineExp = "<< j << "\n";
            for (int b=0; b<splinesExp[j].knots.size()-1; ++b)
                script << nodi[b] << "\n";

            script.close();

            if (allSplineExp == true){

                string pathAndNameString = pathAndPartialName + "knots_p.txt";

                const char* pathAndName = pathAndNameString.c_str();
                ofstream script;
                script.open(pathAndName,std::ios::app);

                script << "knots\t"<<"p = "<< p << "\n";
                for (int b=0; b<splinesExp[p].knots.size()-1; ++b)
                    script << nodi_1[b] << "\n";
                script.close();
            }

            if(allSplineExp == true){ 
                string pathAndNameString = pathAndPartialName + "knots_q.txt";

                const char* pathAndName = pathAndNameString.c_str();
                ofstream script;
                script.open(pathAndName,std::ios::app);

                script << "knots\t"<<"q = "<< q << "\n";
                for (int b=0; b<splinesExp[q].knots.size()-1; ++b)
                    script << nodi_2[b] << "\n";
                script.close();
            }
        }

        if (allSplineExp == true){

            string pathAndNameString = pathAndPartialName + "Exp_D0_p.txt";

            const char* pathAndName = pathAndNameString.c_str();
            ofstream script;
            script.open(pathAndName,std::ios::app);

            script << "x_Exp\t" << fileName << "_Exp\n";
            for (int b=0; b<graphPoints; ++b)
                script << x_1[b] << "\t" << y_D0_1[b] << "\n";

            script.close();
        }

        if (allSplineExp == true){

            string pathAndNameString = pathAndPartialName + "Exp_D0_q.txt";

            const char* pathAndName = pathAndNameString.c_str();
            ofstream script;
            script.open(pathAndName,std::ios::app);

            script << "x_Exp\t" << fileName << "_Exp\n";
            for (int b=0; b<graphPoints; ++b)
                script << x_2[b] << "\t" << y_D0_2[b] << "\n";

            script.close();
        }

    }

}


double Indexes::Stat_SSE(vector <double> b, vector<double> c){

    vector<double> SSE;
    double s;

    for(int i=0; i<b.size(); i++)
        SSE.push_back(pow((b[i]-c[i]),2));

    s = 0;

    for(int i=0; i<SSE.size(); i++)
        s = s + SSE[i];

    return s;   
}