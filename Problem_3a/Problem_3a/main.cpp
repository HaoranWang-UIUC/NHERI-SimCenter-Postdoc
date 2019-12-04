//  Problem_3
//
//  Created by Haoran Wang on 11/27/19.


#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <iomanip>

using namespace std;

void Generate_random_uniform(double* v, int n, double min, double max){
    random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<double> distribution(min,max);
    for (int i=0; i<n; i++)
        v[i] = distribution(generator);
}

void Generate_random_gaussian(double* v, int n, double mean, double std){
    random_device rd;
    default_random_engine generator(rd());
    normal_distribution<double> distribution(mean,std);
    for (int i=0; i<n; i++)
        v[i] = distribution(generator);
}

void Generate_random_parameter(double* rw, double* r, double* Tu, double* Hl, double* Tl, double* Hu, double* L, double* Kw, int num){
    
    Generate_random_uniform(rw, num, 0.1, 0.2);
    Generate_random_gaussian(r, num, 5.0, 2.0);
    Generate_random_gaussian(Tu, num, 0.0, 1.0);//updated later for correlation with Tu
    Generate_random_uniform(Hl, num, 700.0, 820.0);
    Generate_random_gaussian(Tl, num, 0.0, 1.0); //updated later for correlation with Tl
    Generate_random_uniform(L, num, 1120.0, 1680.0);
    Generate_random_uniform(Kw, num, 9500.0, 13000.0);
    double rou = 0.4; //correlation for Tu and Tl
    
    for (int i=0; i<num; i++)
    {
        //Generate Hu conditional on Hl
        r[i] = exp (r[i]);
        double Hu_temp[1];
        Generate_random_uniform(Hu_temp, 1, Hl[i]+100, Hl[i]+200);
        Hu[i] = Hu_temp[0];
        
        //Generate Tu and Tl such that they are correlated
        double Tl_temp = rou * Tu[i] + sqrt(1-rou*rou) * Tl[i];
        Tl[i] = 89.0 + 8.9 * Tl_temp;
        Tu[i] = 89000.0 + 8900.0 * Tu[i];
    }
}

void Calculate_water_flow_rate(double* rw, double* r, double* Tu, double* Hl, double* Tl, double* Hu, double* L, double* Kw, double* y, int num){
    
    const long double pi = 3.141592653589793238L;

    for (int i=0; i<num; i++)
        y[i] = 2 * pi * Tu[i] * (Hu[i] - Hl[i]) / ( log( r[i]/rw[i] ) + 2*L[i]*Tu[i]/rw[i]/rw[i]/Kw[i] + log( r[i]/rw[i] ) * Tu[i]/Tl[i] );

}

//////////////////////////////
//////////////////////////////
int main(int argc, const char * argv[]) {

    int num = 1000;
    double P = 0.0, beta = 150.0;
    double rw[num], r[num], Tu[num], Hl[num], Tl[num], Hu[num], L[num], Kw[num], y[num];
    
    Generate_random_parameter(rw, r, Tu, Hl, Tl, Hu, L, Kw, num);
    
    Calculate_water_flow_rate(rw, r, Tu, Hl, Tl, Hu, L, Kw, y, num);
    
    //Calculate P(y>=beta)
    for (int i=0; i<num; i++)
    {
        if (y[i] >= beta)
            P += 1.0000/num;
    }
    //write results to output file
    ofstream myfile;
    myfile.open ("Random_parameters_and_flow_rate.txt");
    myfile << setw(10) << left << "rw "
           << setw(10) << left << "r "
           << setw(10) << left << "Tu "
           << setw(10) << left << "Hl "
           << setw(10) << left << "Tl "
           << setw(10) << left << "Hu "
           << setw(10) << left << "L "
           << setw(10) << left << "Kw "
           << setw(10) << left<< "y " << endl;
    
    for (int i=0; i<num; i++)
        myfile << setw(10) << left << rw[i]
               << setw(10) << left << r[i]
               << setw(10) << left << Tu[i]
               << setw(10) << left << Hl[i]
               << setw(10) << left << Tl[i]
               << setw(10) << left << Hu[i]
               << setw(10) << left << L[i]
               << setw(10) << left << Kw[i]
               << setw(10) << left << y[i] << endl;
    
    myfile.close();
    cout << "1000 samples have been generated with direct Monte Carlo method!" << endl;
    cout << "Please find results in Random_parameters_and_flow_rate.txt" << endl;
    cout << "The probability for y>=" << beta << " is " << P << endl;

}
