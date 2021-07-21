#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <iomanip>

using namespace std;

#include "Utilities.h"
#include "BasisFunction.h"
#include "Spline.h"
#include "Settings.h"



int indexBestSplineExp;


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

    int index;

    for (int i = 0; i<a.size()-2; i++){
        if (a[i]<=a[i+1] && a[i]<=a[i+2])
            index = i;
        else if (a[i+1]<=a[i] && a[i+1]<=a[i+2])
            index = i+1;
        else if (a[i+2]<=a[i] && a[i+2]<=a[i+1])
            index = i+2;
    }

    return index;
}

// Funziona se splinesExp len is 1?
int calculateBestSpline(vector<double> x, vector<double> y, string criterion, vector<Spline> splinesExp){

    vector<double> ySpl_0;
    vector<double> ySpl_1;
    vector<double> ySpl_2;
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
    
    for (int i=0; i<x.size();i++)
        ySpl_0.push_back(splinesExp[0].D0(x[i]));
    
    for (int i=0; i<x.size();i++)
        ySpl_1.push_back(splinesExp[1].D0(x[i]));

    for (int i=0; i<x.size();i++)
        ySpl_2.push_back(splinesExp[2].D0(x[i]));
    
    SSE.push_back(summedSquaredError(y,ySpl_0));
    SSE.push_back(summedSquaredError(y,ySpl_1));
    SSE.push_back(summedSquaredError(y,ySpl_2));

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

    return indexBestSplineExp;
}

vector<vector<double>> evaluateBestSplineD0 (Spline best_spline) {

    auto x_D0 = vector<double>(graphPoints);
    auto y_D0 = vector<double>(graphPoints);
    vector<vector<double>> splineD0;

    double distance = (best_spline.knots.back()-best_spline.knots[0]) / (double)(graphPoints);

    for (int b=0; b<graphPoints; ++b)
        x_D0[b] = best_spline.knots[0]+(double)b*distance;
    x_D0.back() = best_spline.knots.back();

    // Calculates the ordinates

    for (int b=0; b<graphPoints; ++b)
        y_D0[b] = best_spline.D0(x_D0[b]);

    splineD0.push_back(x_D0);
    splineD0.push_back(y_D0);

    printM(splineD0);

    return splineD0;

}

vector<vector<double>> evaluateBestSplineD1 (Spline best_spline) {

    auto x_D1 = vector<double>(graphPoints);
    auto y_D1 = vector<double>(graphPoints);
    vector<vector<double>> splineD1;

    double distance = (best_spline.knots.back()-best_spline.knots[0]) / (double)(graphPoints);

    for (int b=0; b<graphPoints; ++b)
        x_D1[b] = best_spline.knots[0]+(double)b*distance;
    x_D1.back() = best_spline.knots.back();

    // Calculates the ordinates

    for (int b=0; b<graphPoints; ++b)
        y_D1[b] = best_spline.D1(x_D1[b]);

    splineD1.push_back(x_D1);
    splineD1.push_back(y_D1);

    return splineD1;

}

extern "C"
void edo(double* x, double* y, int length){

    vector<double> x_vector(x, x + length);
    vector<double> y_vector(y, y + length);


    cout <<  "X: ";

    for (int i =0; i < length; i++){
        cout << std::setprecision(12) << x_vector[i] << " ";
    }

    cout << endl;

    cout <<  "Y: ";

    for (int i =0; i < length; i++){
        cout << std::setprecision(12) << y_vector[i] << " ";
    }

    cout << endl;

    vector<Spline> possibleSplines = calculateSplines(x_vector, y_vector, true);

    int index_best = calculateBestSpline(x_vector, y_vector, criterion, possibleSplines);

    cout << "Berst splien: " << indexBestSplineExp << endl;

//    Spline best_spline = possibleSplines[indexBestSplineExp];

    cout << possibleSplines[index_best].D0(2.5) << endl;




    return;
}

int main() {

    vector<vector<double>> splineD0;
    vector<vector<double>> splineD1;

    cout << "Curve Matching\n";

    // Runs the calculations and saves the results

    cout << "\nRunning";

    // qui scrivo io 
    // DATI DI INPUT
    vector<double> x = {0.0 ,0.1,0.2 ,0.3,0.4,0.6,0.8};
    vector<double> y = {0.00021,0.0002,0.00047,0.0005,0.00054,0.00028,0.000205};

    vector<Spline> possibleSplines = calculateSplines(x,y,removeAsymptotes);
    int index_best = calculateBestSpline(x,y,criterion, possibleSplines);
    splineD0 = evaluateBestSplineD0(possibleSplines[index_best]);
    splineD1 = evaluateBestSplineD1(possibleSplines[index_best]);

    printM(splineD0);

    return 0;

}

