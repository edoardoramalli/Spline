class Indexes {

public:

    /* Name of the _exp.txt and _mod.txt files with the data for the splines */
    string fileName;

    ////////////////////////////////////////////////////////////////////////////

    /* Executes all the operations for the comparison of the experimental data
    with the models */
    void solve(const string& folderPath,
               const string& folderName,
               const string& fileName);

    ////////////////////////////////////////////////////////////////////////////

private:

    /* x-coordinates and y-coordinates of the experimental data and the models.
    inputData[0] and inputData[1] contain, respectively, the x-coordinates and
    the y-coordinates of the experimental data, inputData[i] and inputData[i+1]
    contain, respectively, the x-coordinates and the y-coordinates of
    modelNames[i/2-1] */
    vector<vector<double>> inputData;

    /* Splines for the experimental data */
    vector<Spline> splinesExp;

    /* Path to the folder containing the input data */
    string folderPath;

    /* Name of the folder containing the input data */
    string folderName;

    /* Relative errors for each experimental data point */
    vector<double> relativeErrors;

    /* Index in the splinesExp vector of the best fit of the experimental data
    */
    int indexBestSplineExp;

    /* Vector containing the Summed Squared errors of the spline for the experimental data*/
    vector<double> SSE;

    /* Vector containing the likelihood function under the assumption of normal distributed error*/
    vector<double> ll;

    /* Matrix containing AIC AICc BIC for every spline*/
    vector<vector<double>> information;

    /* Degrees of freedom of the spline K in the spline class*/
    vector<int> numOfParam;

    /* Aikake information criterion*/
    vector<double> AIC;
    
    /* Correction of the AIC*/
    vector<double> AICc;

    /* Baesyan information criterion*/
    vector<double> BIC;

    /**/
    vector<double> k; 

    /* Ratio between Dof of the spline and the number of observation used to select between AIC and AICc*/
    vector<double> ratioLK;

    /**/
    vector<double> AICplusAICc;
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

    /* Select the best spline for the experimental data among the three calculated*/
    void calculateIndexBestSplineExp();

    /* Saves the data necessary to plot the splines to .txt files */
    void saveGraphData();
    
    /*Calculate the summed squared error to select the best spline for the 
    experimental data*/

    double summedSquaredError(vector <double> b, vector<double> c);

    /* Function to calculate the likelihood function under the hypothesis of normal distributed
        errors*/
    vector<double> logLikeliHood(double n, vector<double> residuals);

    /* Function to calculate AIC AICc BIC*/
    vector<vector<double>> informationCriterion(vector<double> ll, double n, vector<int> numOfParam, vector<double> residuals);

    /* Function to find the minimum of a vector of size 3*/
    int positionOfMinimum(vector<double> a);

};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



void Indexes::solve(const string& FolderPath,const string& FolderName,const string& FileName) {

    folderPath = FolderPath;
    folderName = FolderName;
    fileName = FileName;

    // Processes the data in the input files
    this->readData();

    // Calculates the splines for the experimental data
    this->calculateSplines();

    // Calculates the scores from the indexes
    this->calculateIndexBestSplineExp();

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
    if (removeAsymptotes == true){
        splinesExp[0].removeAsymptotes();
    }
    splinesExp[0].normalizeCoefficients(
                                    -splinesExp[0].yD0Min,
                                    splinesExp[0].yD0Max-splinesExp[0].yD0Min,
                                    splinesExp[0].yD1MaxAbs);


    if (splinesExp.size() > 1) {
        splinesExp[1].solve(inputData[0],inputData[1],0,2);
        if (removeAsymptotes == true){
            splinesExp[1].removeAsymptotes();
        }
        splinesExp[1].normalizeCoefficients(
                                    -splinesExp[1].yD0Min,
                                    splinesExp[1].yD0Max-splinesExp[1].yD0Min,
                                    splinesExp[1].yD1MaxAbs);

    }

    if (splinesExp.size() > 2) {
        splinesExp[2].solve(inputData[0],inputData[1],0,5);
        if (removeAsymptotes == true){
            splinesExp[2].removeAsymptotes();
        }
        splinesExp[2].normalizeCoefficients(
                                    -splinesExp[2].yD0Min,
                                    splinesExp[2].yD0Max-splinesExp[2].yD0Min,
                                    splinesExp[2].yD1MaxAbs);

    }

}

void Indexes::calculateIndexBestSplineExp() {

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
    
    SSE.push_back(summedSquaredError(inputData[1],ySpl_0));
    SSE.push_back(summedSquaredError(inputData[1],ySpl_1));
    SSE.push_back(summedSquaredError(inputData[1],ySpl_2));

    double numOfObs = inputData[0].size();

    ll = logLikeliHood(numOfObs,SSE);

    for(int i=0;i<splinesExp.size();i++)
        numOfParam.push_back(splinesExp[i].K);

    information = informationCriterion(ll,numOfObs,numOfParam,SSE);

    AIC = information[0];
    AICc = information[1];
    BIC = information [2];
    k = information [3];

    for (int i=0; i<k.size();i++)
        ratioLK.push_back(k[i]/numOfObs);

    if (criterion == "SSE"){

        indexBestSplineExp = positionOfMinimum(SSE);
    
    }

    if (criterion == "AIC"){

    	for (int i=0; i<ratioLK.size(); i++){
    		if (ratioLK[i] <= numberOfRatiolkForAICcUse)
    			AICplusAICc.push_back(AICc[i]);
    		else
    			AICplusAICc.push_back(AIC[i]);
    	}   

        indexBestSplineExp = positionOfMinimum(AICplusAICc);
    }

    if (criterion == "BIC"){
        
        indexBestSplineExp = positionOfMinimum(BIC);
    }

}


void Indexes::saveGraphData() {

    if (graphs == false) return;

    int j = indexBestSplineExp;
    int p;
    int q;

    auto x = vector<double>(graphPoints);
    auto x_1 = vector<double>(graphPoints);
    auto x_2 = vector<double>(graphPoints);

    vector<double> y_D0, y_D1, nodi;
    vector<double> y_D0_1, nodi_1;
    vector<double> y_D0_2, nodi_2;
    if (graphsD0 == true) y_D0 = vector<double>(graphPoints);
    if (graphsD1 == true) y_D1 = vector<double>(graphPoints);
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

        if (coseUtili == true){

            string pathAndNameString = pathAndPartialName + "recap.txt";

            const char* pathAndName = pathAndNameString.c_str();
            ofstream script;
            script.open(pathAndName,std::ios::app);

            script <<  "DoF"<< setw(10) << "index"<< setw(10) << "SSE" << setw(15)<< "AIC" << setw(14)<< "AICc" << setw(10)<< "BIC" << setw(15) << "ratioLK" << setw(15)<< "AICorAICc\n";
            script << splinesExp[j].K << setw(7) << j << setw(22) << SSE[j]<< setw(12) << AIC[j] << setw(12) << AICc[j] << setw(10) << BIC[j]<<setw(15)<<ratioLK[j]<<setw(20)<<AICplusAICc[j]<<"\n";
            script << splinesExp[p].K << setw(7) << p << setw(22) << SSE[p]<< setw(12) << AIC[p] << setw(12) << AICc[p] << setw(10) << BIC[p]<<setw(15)<<ratioLK[p]<<setw(20)<<AICplusAICc[p]<<"\n";
            script << splinesExp[q].K << setw(7) << q << setw(22) << SSE[q]<< setw(12) << AIC[q] << setw(12) << AICc[q] << setw(10) << BIC[q]<<setw(15)<<ratioLK[q]<<setw(20)<<AICplusAICc[q]<<"\n";

            script.close();
        }

    }

}

double Indexes::summedSquaredError(vector <double> b, vector<double> c){

    vector<double> SSE;
    double s;

    for(int i=0; i<b.size(); i++)
        SSE.push_back(pow((b[i]-c[i]),2));

    s = 0;

    for(int i=0; i<SSE.size(); i++)
        s = s + SSE[i];

    return s;  

}

vector<double> Indexes::logLikeliHood(double n, vector<double> residuals){

    vector<double> ll;

    for(int i=0; i<SSE.size(); i++)
        ll.push_back(n * log(residuals[i] / n));

    return ll; 

}

vector<vector<double>> Indexes::informationCriterion(vector<double> ll, double n, vector<int> numOfParam,
                                                         vector<double> residuals){
    
    vector<vector<double>> information ;
    vector<double> AIC;
    vector<double> correctionAIC;
    vector<double> AICc;
    vector<double> BIC;
    vector<double> k; 

    for(int i=0;i<numOfParam.size();i++)
        k.push_back(2*(numOfParam[i]+1)+1);

    for(int i=0; i<ll.size();i++){
        
        AIC.push_back(2*k[i]+ll[i]);
        correctionAIC.push_back(2*k[i]*(k[i]+1)/(n-k[i]-1));
        BIC.push_back(ll[i]+k[i]*log(n));
    }
    for (int i = 0; i<correctionAIC.size(); i++)
        AICc.push_back(AIC[i]+correctionAIC[i]);

    information.push_back(AIC);
    information.push_back(AICc);
    information.push_back(BIC);
    information.push_back(k);
    
    return information;
        
}

int Indexes::positionOfMinimum(vector<double> a){

    int index;

    for (int i = 0; i<AIC.size()-2; i++){
        if (a[i]<=a[i+1] && a[i]<=a[i+2])
            index = i;
        else if (a[i+1]<=a[i] && a[i+1]<=a[i+2])
            index = i+1;
        else if (a[i+2]<=a[i] && a[i+2]<=a[i+1])
            index = i+2;
    }

    return index;
}