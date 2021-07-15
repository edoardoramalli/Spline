/* Reads the program settings from Setting.txt */
void readProgramSettings() {

    string possibleNegativeOrdinatesString;
    string lineUpMaximaString = "yes";
    string useSumOfIndexesForAlignmentString = "yes";
    string saveIndexesToCsvString;
    string graphsFittingD0String;
    string graphsD0String;
    string graphsD1String;
    string graphsD0ShiftString;
    string graphsD1ShiftString;
    // aggiunte da Timoteo Dinelli per confrontare le diverse spline generate e salvare 
    // i dati per la generazione dei grafi in un modo differente
    string allSplineExpString;
    string knotsToSaveString;
    string saveExpDataString;

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
        else if (s == "@fractionOfExpHeightForZeroLineIdentification:")
            iss >> fractionOfExpHeightForZeroLineIdentification;
        else if (s == "@fractionOfExpRangeForModelExtrapolation:")
            iss >> fractionOfExpRangeForModelExtrapolation;
        else if (s == "@fractionOfExpRangeForMinShift:")
            iss >> fractionOfExpRangeForMinShift;
        else if (s == "@fractionOfExpRangeForShiftAroundMaximum:")
            iss >> fractionOfExpRangeForShiftAroundMaximum;
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
        else if (s == "@saveIndexesToCsv:")
            iss >> saveIndexesToCsvString;
        else if (s == "@graphsFittingD0:")
            iss >> graphsFittingD0String;
        else if (s == "@graphsD0:")
            iss >> graphsD0String;
        else if (s == "@graphsD1:")
            iss >> graphsD1String;
        else if (s == "@graphsD0Shift:")
            iss >> graphsD0ShiftString;
        else if (s == "@graphsD1Shift:")
            iss >> graphsD1ShiftString;
        else if (s == "@graphPoints:")
            iss >> graphPoints;
        else if (s == "@allSplineExp:")
            iss >> allSplineExpString;
        else if (s == "@knotsToSave:")
            iss >> knotsToSaveString;
        else if (s == "@saveExpData:")
            iss >> saveExpDataString;

    }
    myfile.close();

    if (possibleNegativeOrdinatesString == "yes")
        possibleNegativeOrdinates = true;

    if (lineUpMaximaString != "yes")
        lineUpMaxima = false;

    if (useSumOfIndexesForAlignmentString != "yes")
        useSumOfIndexesForAlignment = false;

    if (saveIndexesToCsvString == "yes")
        saveIndexesToCsv = true;

    if (graphsFittingD0String != "yes") graphsFittingD0 = false;
    if (graphsD0String != "yes") graphsD0 = false;
    if (graphsD1String == "yes") graphsD1 = true;
    if (graphsD0ShiftString != "yes") graphsD0Shift = false;
    if (graphsD1ShiftString == "yes") graphsD1Shift = true;
    if (allSplineExpString == "yes") allSplineExp = true;
    if (knotsToSaveString == "yes") knotsToSave = true;
    if (saveExpDataString == "yes") saveExpData = true;

    if (graphsFittingD0 == false)
        if (graphsD0 == false)
            if (graphsD1 == false)
                if (graphsD0Shift == false)
                    if (graphsD1Shift == false)
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
    if (fractionOfExpHeightForZeroLineIdentification != 0.02)
        s += "fractionOfExpHeightForZeroLineIdentification: "+
             to_string(fractionOfExpHeightForZeroLineIdentification)+"\n";
    if (fractionOfExpRangeForModelExtrapolation != 0.5)
        s += "fractionOfExpRangeForModelExtrapolation: "+
             to_string(fractionOfExpRangeForModelExtrapolation)+"\n";
    if (fractionOfExpRangeForMinShift != 0.005)
        s += "fractionOfExpRangeForMinShift: "+
             to_string(fractionOfExpRangeForMinShift)+"\n";
    if (fractionOfExpRangeForShiftAroundMaximum != 0.05)
        s += "fractionOfExpRangeForShiftAroundMaximum: "+
             to_string(fractionOfExpRangeForShiftAroundMaximum)+"\n";
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
    if (saveIndexesToCsv != false)
        s += "saveIndexesToCsv: true\n";
    if (knotsToSave == true)
        s += "Saving knots for spline calculations\n";
    if (allSplineExp == true)
        s += "Saving the data for every spline for the experimental data \n";

    if (s != "")
        cout << "\nNon-default settings:\n" << s;


    

    if (graphs == true) {

        string graphString = "\nGraphs: ";
        if (graphsFittingD0 == true)
            graphString += "FittingD0";
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
        if (graphsD0Shift == true && useSumOfIndexesForAlignment == true) {
            if (graphString != "\nGraphs: ")
                graphString += ", ";
            graphString += "D0Shift";
        }
        if (graphsD1Shift == true && useSumOfIndexesForAlignment == true) {
            if (graphString != "\nGraphs: ")
                graphString += ", ";
            graphString += "D1Shift";
        }
        cout << graphString << "\n";

        string graphSettings = "";

        if (graphPoints != 500)
            graphSettings += "graphPoints: " + to_string(graphPoints);

        if (graphSettings != "")
            cout << graphSettings << "\n";

    }

}