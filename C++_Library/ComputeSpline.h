vector<Spline> calculateSplines(vector<double> x, vector<double> y, bool removeAsymptotes) {

    vector<Spline> splinesExp;

    if (x.size() < 3)
        splinesExp = vector<Spline>(1);
    else if (x.size() < 5)
        splinesExp = vector<Spline>(2);
    else
        splinesExp = vector<Spline>(3);

    splinesExp[0].solve(x,y,0,0);
    if (removeAsymptotes == true){
        splinesExp[0].removeAsymptotes();
    }
    splinesExp[0].normalizeCoefficients(
                                    -splinesExp[0].yD0Min,
                                    splinesExp[0].yD0Max-splinesExp[0].yD0Min,
                                    splinesExp[0].yD1MaxAbs);


    if (splinesExp.size() > 1) {
        splinesExp[1].solve(x,y,0,2);
        if (removeAsymptotes == true){
            splinesExp[1].removeAsymptotes();
        }
        splinesExp[1].normalizeCoefficients(
                                    -splinesExp[1].yD0Min,
                                    splinesExp[1].yD0Max-splinesExp[1].yD0Min,
                                    splinesExp[1].yD1MaxAbs);

    }

    if (splinesExp.size() > 2) {
        splinesExp[2].solve(x,y,0,5);
        if (removeAsymptotes == true){
            splinesExp[2].removeAsymptotes();
        }
        splinesExp[2].normalizeCoefficients(
                                    -splinesExp[2].yD0Min,
                                    splinesExp[2].yD0Max-splinesExp[2].yD0Min,
                                    splinesExp[2].yD1MaxAbs);

    }

    return splinesExp;

}

double summedSquaredError(vector <double> b, vector<double> c){

    vector<double> SSE;
    double s;

    for(int i=0; i<b.size(); i++)
        SSE.push_back(pow((b[i]-c[i]),2));

    s = 0;

    for(int i=0; i<SSE.size(); i++)
        s = s + SSE[i];

    return s;

}

vector<double> logLikeliHood(double n, vector<double> residuals){

    vector<double> ll;

    for(int i=0; i<residuals.size(); i++)
        ll.push_back(n * log(residuals[i] / n));

    return ll;

}

vector<vector<double>> informationCriterion(vector<double> ll, double n, vector<int> numOfParam,vector<double> residuals){

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

int positionOfMinimum(vector<double> a){

    int indexmin = 0;

    for (unsigned i = 0; i < a.size(); ++i)
    {
        if (a[i] <= a[indexmin]) // Found a smaller min
            indexmin = i;
    }

    return indexmin;
}

int calculateBestSpline(vector<double> x, vector<double> y, string criterion, vector<Spline> splinesExp){

//    vector<double> ySpl_0;
//    vector<double> ySpl_1;
//    vector<double> ySpl_2;
//

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
    double numOfObs = x.size();

    int indexBestSplineExp;

    for (int k=0; k < splinesExp.size(); k++){
        vector<double> ySpl_tmp;
        for (int i=0; i<x.size();i++)
            ySpl_tmp.push_back(splinesExp[k].D0(x[i]));
        SSE.push_back(summedSquaredError(y,ySpl_tmp));
    }

//    if (splinesExp.size()==1){
//
//        for (int i=0; i<x.size();i++)
//            ySpl_0.push_back(splinesExp[0].D0(x[i]));
//        SSE.push_back(summedSquaredError(y,ySpl_0));
//    }
//    else if (splinesExp.size()==2){
//
//        for (int i=0; i<x.size();i++){
//            ySpl_0.push_back(splinesExp[0].D0(x[i]));
//        }
//
//
//        for (int i=0; i<x.size();i++){
//            ySpl_1.push_back(splinesExp[1].D0(x[i]));
//        }
//
//
//        SSE.push_back(summedSquaredError(y,ySpl_0));
//
//        SSE.push_back(summedSquaredError(y,ySpl_1));
//    }
//    else{
//
//        for (int i=0; i<x.size();i++)
//            ySpl_0.push_back(splinesExp[0].D0(x[i]));
//
//        for (int i=0; i<x.size();i++)
//            ySpl_1.push_back(splinesExp[1].D0(x[i]));
//
//        for (int i=0; i<x.size();i++)
//            ySpl_2.push_back(splinesExp[2].D0(x[i]));
//
//        SSE.push_back(summedSquaredError(y,ySpl_0));
//        SSE.push_back(summedSquaredError(y,ySpl_1));
//        SSE.push_back(summedSquaredError(y,ySpl_2));
//    }

    ll = logLikeliHood(numOfObs,SSE);

    for(int i=0;i<splinesExp.size();i++){
        numOfParam.push_back(splinesExp[i].K);
    }


    information = informationCriterion(ll,numOfObs,numOfParam,SSE);

    AIC = information[0];
    AICc = information[1];
    BIC = information[2];
    k = information[3];

    for (int i=0; i<k.size();i++){
        ratioLK.push_back(k[i]/numOfObs);
    }


    if (criterion == "SSE"){
        indexBestSplineExp = positionOfMinimum(SSE);
    }

    if (criterion == "AIC"){
        for (int i=0; i<ratioLK.size(); i++){
            if (ratioLK[i] <= numberOfRatiolkForAICcUse){
                AICplusAICc.push_back(AICc[i]);
            }
            else{
                AICplusAICc.push_back(AIC[i]);
            }
        }
        indexBestSplineExp = positionOfMinimum(AICplusAICc);
    }

    if (criterion == "BIC"){
        indexBestSplineExp = positionOfMinimum(BIC);
    }

    return indexBestSplineExp;
}
