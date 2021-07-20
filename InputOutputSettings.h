/* Reads the program settings from Setting.txt */
void readProgramSettings() {

    string possibleNegativeOrdinatesString;
    string graphsD0String;
    string graphsD1String;
    // aggiunte da Timoteo Dinelli per confrontare le diverse spline generate e salvare 
    // i dati per la generazione dei grafi in un modo differente
    string allSplineExpString;
    string knotsToSaveString;
    string saveExpDataString;
    string coseUtiliString;
    string criterionString;
    string removeAsymptotesString;

    string line;
    ifstream myfile("./Settings.txt");
    while (getline(myfile,line)) {

        istringstream iss(line);
        string s;
        iss >> s;

        if (s == "@lambdaSearchInterval:")
            iss >> lambdaSearchInterval;
        else if (s == "@numberOfStepsLambda:")
            iss >> numberOfStepsLambda;
        else if (s == "@fractionOfOrdinateRangeForAsymptoteIdentification:")
            iss >> fractionOfOrdinateRangeForAsymptoteIdentification;
        else if (s == "@fractionOfOrdinateRangeForMaximumIdentification:")
            iss >> fractionOfOrdinateRangeForMaximumIdentification;
        else if (s == "@possibleNegativeOrdinates:")
            iss >> possibleNegativeOrdinatesString;
        else if (s == "@numberOfTrapezoids:")
            iss >> numberOfTrapezoids;
        else if (s == "@defaultRelativeError:")
            iss >> defaultRelativeError;
        else if (s == "@graphsD0:")
            iss >> graphsD0String;
        else if (s == "@graphsD1:")
            iss >> graphsD1String;
        else if (s == "@graphPoints:")
            iss >> graphPoints;
        else if (s == "@allSplineExp:")
            iss >> allSplineExpString;
        else if (s == "@knotsToSave:")
            iss >> knotsToSaveString;
        else if (s == "@saveExpData:")
            iss >> saveExpDataString;
        else if (s == "@coseUtili:")
        	iss >> coseUtiliString;
        else if (s == "@criterion:")
            iss >> criterionString;
        else if (s == "@removeAsymptotes:")
            iss >> removeAsymptotesString;
        else if (s == "@numberOfRatiolkForAICcUse:")
            iss >> numberOfRatiolkForAICcUse;
    }
    myfile.close();

    if (possibleNegativeOrdinatesString == "yes")
        possibleNegativeOrdinates = true;

    if (graphsD0String != "yes") graphsD0 = false;
    if (graphsD1String == "yes") graphsD1 = true;
    if (allSplineExpString == "yes") allSplineExp = true;
    if (knotsToSaveString == "yes") knotsToSave = true;
    if (saveExpDataString == "yes") saveExpData = true;
    if (coseUtiliString == "yes") coseUtili = true;
    if (criterionString == "AIC") criterion = "AIC";
    if (criterionString == "SSE") criterion = "SSE";
    if (criterionString == "BIC") criterion = "BIC";
    if (removeAsymptotesString == "no") removeAsymptotes = false;

    if (graphsD0 == false)
        if (graphsD1 == false)
            graphs = false;

}



/* Displays the non-default Curve Matching settings and the graph settings */
void displaySettings() {

    string s = "";

    if (lambdaSearchInterval != 6)
        s += "lambdaSearchInterval: "+to_string(lambdaSearchInterval)+"\n";
    if (numberOfStepsLambda != 13)
        s += "numberOfStepsLambda: "+to_string(numberOfStepsLambda)+"\n";
    if (numberOfRatiolkForAICcUse != 40)
        s += "numberOfRatiolkForAICcUse: "+to_string(numberOfRatiolkForAICcUse)+"\n";
    if (fractionOfOrdinateRangeForAsymptoteIdentification != 0.005)
        s += "fractionOfOrdinateRangeForAsymptoteIdentification: "+
             to_string(fractionOfOrdinateRangeForAsymptoteIdentification)+"\n";
    if (fractionOfOrdinateRangeForMaximumIdentification != 0.025)
        s += "fractionOfOrdinateRangeForMaximumIdentification: "+
             to_string(fractionOfOrdinateRangeForMaximumIdentification)+"\n";
    if (possibleNegativeOrdinates != false)
        s += "possibleNegativeOrdinates: true\n";
    if (numberOfTrapezoids != 99)
        s += "numberOfTrapezoids: "+to_string(numberOfTrapezoids)+"\n";
    if (defaultRelativeError != 0.1)
        s += "defaultRelativeError: "+to_string(defaultRelativeError)+"\n";
    if (knotsToSave == true)
        s += "Saving knots for spline calculations\n";
    if (allSplineExp == true)
        s += "Saving the data for every spline for the experimental data \n";
    if (coseUtili == true)
        s += "Generation of the recap file \n";
    if (criterion == "SSE")
        s += "The choice of the best spline for the experimental data is done minimizing the SSE \n";
    if (criterion == "BIC")
        s += "The choice of the best spline for the experimental data is done minimizing the BIC \n";
    if (removeAsymptotes == false)
        s += "Not removing the left and right asymptotes for the spline";

    if (s != "")
        cout << "\nNon-default settings:\n" << s;


    

    if (graphs == true) {

        string graphString = "\nGraphs: ";
        if (graphsD0 == true) {
            if (graphString != "\nGraphs: ")
                graphString += ", ";
            graphString += "D0";
        }
        if (graphsD1 == true) {
            if (graphString != "\nGraphs: ")
                graphString += ", ";
            graphString += "D1";
        }
        
        cout << graphString << "\n";

        string graphSettings = "";

        if (graphPoints != 500)
            graphSettings += "graphPoints: " + to_string(graphPoints);

        if (graphSettings != "")
            cout << graphSettings << "\n";

    }

}