/* Generates the splines according to the splineTypes and number of points */
vector<Spline> calculateSplines(vector<double> x, vector<double> y, int splineType) {

    vector<Spline> splines(3);

    vector<int> numberOfAbscissaeSeparatingConsecutiveKnots_vector = {0, 2, 5};

    if (splineType == 1 || x.size() < 3){  // if model or len(x) < 3 --> just one spline
        splines.resize(1);
    }
    else if (x.size() < 5){
        splines.resize(2);
    }
    else{
        splines.resize(3);
    }

    for (int i = 0; i < (int)splines.size(); i++) {
        splines[i].solve(x, y, splineType, numberOfAbscissaeSeparatingConsecutiveKnots_vector[i]);
    }

    return splines;

}

double summedSquaredError(vector <double> b, vector<double> c){
    double SSE = 0;
    for(int i=0; i < (int)b.size(); i++){
        SSE += pow((b[i] - c[i]), 2);
    }
    return SSE;
}

vector<double> logLikeliHood(double n, vector<double> residuals){

    vector<double> ll;

    for(int i=0; i < (int)residuals.size(); i++)
        ll.push_back(n * log(residuals[i] / n));

    return ll;

}

vector<vector<double>> informationCriterion(vector<double> ll, int n, vector<int> numOfParam){

    vector<vector<double>> information ;
    vector<double> AIC;
    vector<double> correctionAIC;
    vector<double> AICc;
    vector<double> BIC;
    vector<double> k;

    for(int i=0; i < (int)numOfParam.size();i++)
        k.push_back(2*(numOfParam[i]+1)+1);

    for(int i=0; i < (int)ll.size();i++){

        AIC.push_back(2*k[i]+ll[i]);
        correctionAIC.push_back(2*k[i]*(k[i]+1)/(n-k[i]-1));
        BIC.push_back(ll[i]+k[i]*log(n));
    }
    for (int i = 0; i < (int)correctionAIC.size(); i++)
        AICc.push_back(AIC[i]+correctionAIC[i]);

    information.push_back(AIC);
    information.push_back(AICc);
    information.push_back(BIC);
    information.push_back(k);

    return information;

}

/* Return the index of the minimum element */
int positionOfMinimum(vector<double> v){
    return min_element(v.begin(), v.end()) - v.begin();
}

/* Given a vector of Splines return the best spline based on the criterion */
int calculateBestSpline(vector<Spline> splines, string criterion){

    // If the length of splines is 1 then the only spline is the best spline
    if (splines.size() == 1){
        return 0;
    }

    vector<int> numOfParam;
    vector<double> AIC;
    vector<double> AICc;
    vector<double> BIC;
    vector<double> AICplusAICc;
    vector<double> SSE;
    vector<double> ll;
    vector<vector<double>> information;
    vector<double> ratioLK;
    vector<double> k;
    // The original X vector is in every Spline. I take it from the first one
    int numOfObs = splines[0].originalAbscissae.size();

    int indexBestSpline;

    for (int k=0; k < (int)splines.size(); k++){
        vector<double> ySpl_tmp;
        for (int i=0; i < numOfObs;i++){
            ySpl_tmp.push_back(splines[k].D0(splines[0].originalAbscissae[i]));
        }
        SSE.push_back(summedSquaredError(splines[0].originalOrdinates, ySpl_tmp));
    }

    ll = logLikeliHood(numOfObs,SSE);

    for(int i=0; i< (int)splines.size(); i++){
        numOfParam.push_back(splines[i].K);
    }

    information = informationCriterion(ll, numOfObs, numOfParam);

    AIC = information[0];
    AICc = information[1];
    BIC = information[2];
    k = information[3];


    if (criterion == "SSE"){
        indexBestSpline = positionOfMinimum(SSE);
    }

    if (criterion == "AIC"){
        for (int i=0; i < (int)k.size();i++){
            ratioLK.push_back(k[i] / numOfObs);
        }
        for (int i=0; i < (int)ratioLK.size(); i++){
            if (ratioLK[i] <= numberOfRatiolkForAICcUse){
                AICplusAICc.push_back(AICc[i]);
            }
            else{
                AICplusAICc.push_back(AIC[i]);
            }
        }
        indexBestSpline = positionOfMinimum(AICplusAICc);
    }


    if (criterion == "BIC"){
        indexBestSpline = positionOfMinimum(BIC);
    }

    return indexBestSpline;
}

vector<vector<double>> evaluateSpline (Spline best_spline, int der) {

    auto x_eval = vector<double>(graphPoints);
    auto y_eval = vector<double>(graphPoints);
    vector<vector<double>> spline_evaluate;


    double (Spline::*evaluate_function)(double);

    if (der == 0){
        evaluate_function = &Spline::D0;
    }
    else if (der == 1){
        evaluate_function = &Spline::D1;
    }
    else if (der == 2){
        evaluate_function = &Spline::D2;
    }

    double distance = (best_spline.knots.back()-best_spline.knots[0]) / (double)(graphPoints);

    for (int b=0; b<graphPoints; ++b){
        x_eval[b] = best_spline.knots[0]+(double)b*distance;
    }

    x_eval.back() = best_spline.knots.back();

    // Calculates the ordinates

    for (int b=0; b<graphPoints; ++b){
        y_eval[b] = (best_spline.*evaluate_function)(x_eval[b]);
    }


    spline_evaluate.push_back(x_eval);
    spline_evaluate.push_back(y_eval);

    return spline_evaluate;

}
