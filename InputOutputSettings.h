/* Reads the program settings from Setting.txt */
void readProgramSettings() {

    string possibleNegativeOrdinatesString;
    string lineUpMaximaString = "yes";
    string useSumOfIndexesForAlignmentString = "yes";
    string graphsD0String;
    string graphsD1String;
    // aggiunte da Timoteo Dinelli per confrontare le diverse spline generate e salvare 
    // i dati per la generazione dei grafi in un modo differente
    string allSplineExpString;
    string knotsToSaveString;
    string saveExpDataString;
    string coseUtiliString;
    string criterionSSEString;
    string criterionAICString;
    string criterionBICString;
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
        else if (s == "@numberOfBootstrapVariations:")
            iss >> numberOfBootstrapVariations;
        else if (s == "@lineUpMaxima:")
            iss >> lineUpMaximaString;
        else if (s == "@useSumOfIndexesForAlignment:")
            iss >> useSumOfIndexesForAlignmentString;
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
        else if (s == "@criterionSSE:")
            iss >> criterionSSEString;
        else if (s == "@criterionAIC:")
            iss >> criterionAICString;
        else if (s == "@criterionBIC:")
            iss >> criterionBICString;
        else if (s == "@removeAsymptotes:")
            iss >> removeAsymptotesString;

    }
    myfile.close();

    if (possibleNegativeOrdinatesString == "yes")
        possibleNegativeOrdinates = true;

    if (lineUpMaximaString != "yes")
        lineUpMaxima = false;

    if (useSumOfIndexesForAlignmentString != "yes")
        useSumOfIndexesForAlignment = false;

    if (graphsD0String != "yes") graphsD0 = false;
    if (graphsD1String == "yes") graphsD1 = true;
    if (allSplineExpString == "yes") allSplineExp = true;
    if (knotsToSaveString == "yes") knotsToSave = true;
    if (saveExpDataString == "yes") saveExpData = true;
    if (coseUtiliString == "yes") coseUtili = true;
    if (criterionSSEString == "yes") criterionSSE = true, criterionAIC = false, criterionBIC = false;
    if (criterionAICString == "yes") criterionSSE = false, criterionAIC = true, criterionBIC = false;
    if (criterionBICString == "yes") criterionSSE = false, criterionAIC = false, criterionBIC = true;
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
    if (fractionOfOrdinateRangeForAsymptoteIdentification != 0.005)
        s += "fractionOfOrdinateRangeForAsymptoteIdentification: "+
             to_string(fractionOfOrdinateRangeForAsymptoteIdentification)+"\n";
    if (fractionOfOrdinateRangeForMaximumIdentification != 0.025)
        s += "fractionOfOrdinateRangeForMaximumIdentification: "+
             to_string(fractionOfOrdinateRangeForMaximumIdentification)+"\n";
    if (possibleNegativeOrdinates != false)
        s += "possibleNegativeOrdinates: true\n";
    if (numberOfBootstrapVariations != 20)
        s += "numberOfBootstrapVariations: "+
             to_string(numberOfBootstrapVariations)+"\n";
    if (lineUpMaxima != true)
        s += "lineUpMaxima: false\n";
    if (useSumOfIndexesForAlignment != true)
        s += "useSumOfIndexesForAlignment: false\n";
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
    if (criterionSSE == true)
        s += "The choice of the best spline for the experimental data is done minimizing the SSE \n";
    if (criterionBIC == true)
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