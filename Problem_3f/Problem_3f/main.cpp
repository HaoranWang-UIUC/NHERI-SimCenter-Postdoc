//  Problem_3f
//  Created by Haoran Wang on 12/4/19.

#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
//#include <iomanip>
#include <algorithm>
//#include <sstream>
//#include <map>

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

void Generate_random_parameter(double* rw, double* logr, double* Tu, double* Tl, double* L, double* Kw, double* sigma_e, int num){
    
    Generate_random_uniform(rw, num, 0.1, 0.2);
    Generate_random_gaussian(logr, num, 5.0, 2.0);
    Generate_random_gaussian(Tu, num, 0.0, 1.0);//updated later for correlation with Tu
    Generate_random_gaussian(Tl, num, 0.0, 1.0); //updated later for correlation with Tl
    Generate_random_uniform(L, num, 1120.0, 1680.0);
    Generate_random_uniform(Kw, num, 9500.0, 13000.0);
    Generate_random_uniform(sigma_e, num, 0.0, 10.0);
    double rou = 0.4; //correlation for Tu and Tl
    
    for (int i=0; i<num; i++)
    {
        //Generate Tu and Tl such that they are correlated by rou
        double Tl_temp = rou * Tu[i] + sqrt(1-rou*rou) * Tl[i];
        Tl[i] = 89.0 + 8.9 * Tl_temp;
        Tu[i] = 89000.0 + 8900.0 * Tu[i];
    }
}

void Read_input(double* u_m, double* y_m)
{
    fstream file1, file2;
    file1.open("input.txt", ios::in);
    file2.open("output.txt", ios::in);
    int i = 0;
    while (file1.good())
    {
        file1 >> u_m[i];
        i++;
    }
    
    i = 0;
    while (file2.good())
    {
        file2 >> y_m[i];
        i++;
    }
}

double Flow_rate(double u, double rw, double logr, double Tu, double Tl, double L, double Kw)
{
    const long double pi = 3.141592653589793238L;
    double y;
    y = 2 * pi * Tu * u / ( (logr-log(rw)) + 2*L*Tu/rw/rw/Kw + ( logr-log(rw) ) * Tu/Tl );
    return y;
}

long double likely_function(double sigma_e, int M, double error)
{
    const long double pi = 3.141592653589793238L;
    long double p;
    p = exp(-1.0 * error/2.0/sigma_e/sigma_e);
    p = p * 1.0e40 / pow(2.0*pi*sigma_e*sigma_e, 20);
    return p;
}

double max_posterior(double* p, int num)
{
    for (int i=0; i<num; i++)
        if (isnan(p[i])) p[i]=0.0;
    
    return *max_element(p, p+num);
}

//Rejection method for sampling on a pdf
void Generate_sample(double* sample, int num_sample, double* p, double min, double max)
{
    for (int i=0; i<num_sample; i++)
    {
        int flag=0;
        double step = (max-min)/1000.0, p_max=max_posterior(p,1001);
        while (flag==0)
        {
            double temp[1], p_temp=0, ruler[1];
            Generate_random_uniform(temp, 1, min, max); //to speed up the sampling, I decrease the domain based on its pdf
            Generate_random_uniform(ruler, 1, 0.0, 1.0);
            int j = floor ( (temp[0] - min)/step );
            p_temp =(p[j+1]-p[j])/step * (temp[0] - min-j*step) + p[j];
            
            if (p_temp/p_max >= ruler[0])
            {   flag=1;
                sample[i] = temp[0]; }
        }
    }
}
void Sample_p_sigma_e(double* sigma_e_new, double* rw, double* logr, double* Tu, double* Tl, double* L, double* Kw, double* u_m, double* y_m, int M, int num_mc, int num_sample)
{
    double sigma_e, error[num_mc], error_temp;
    for (int j=0; j<num_mc; j++)
        for (int i=0; i<M; i++)
        {
            error_temp = y_m[i] - Flow_rate(u_m[i], rw[j], logr[j], Tu[j], Tl[j], L[j], Kw[j]);
            error[j] += error_temp * error_temp;
        }
    double p_sigma_e[1001];
    for (int i=0; i<1001; i++)
    {   p_sigma_e[i] = 0;
        for (int j=0; j<num_mc; j++)
        {
            sigma_e = 0.0 + 0.01*i;
            p_sigma_e[i] += likely_function(sigma_e, M, error[j]); //calculate posterior sitribution
        }
    }
    
    Generate_sample(sigma_e_new, num_sample, p_sigma_e, 1.0, 5.0);//Rejection method for sampling on a pdf
    cout << "New samples for sigma_e have been generated!" << endl;
}

void Sample_p_rw(double* sigma_e, double* rw_new, double* logr, double* Tu, double* Tl, double* L, double* Kw, double* u_m, double* y_m, int M, int num_mc, int num_sample)
{
    double rw, error[num_mc], error_temp, p_rw[1001];
    
    for (int ii=0; ii<1001; ii++)
    {   rw = 0.1 + ii*0.0001;
        p_rw[ii] = 0;
        
        for (int j=0; j<num_mc; j++)
        {   error[j] = 0.0;
            
            for (int i=0; i<M; i++)
            {
                error_temp = y_m[i] - Flow_rate(u_m[i], rw, logr[j], Tu[j], Tl[j], L[j], Kw[j]);
                error[j] += error_temp * error_temp;
            }
            
            p_rw[ii] += likely_function(sigma_e[j], M, error[j]);
        }
    }

       Generate_sample(rw_new, num_sample, p_rw, 0.1, 0.2);//Rejection method for sampling on a pdf
       cout << "New samples for rw have been generated!" << endl;
}

void Sample_p_r(double* sigma_e, double* rw, double* r_new, double* Tu, double* Tl, double* L, double* Kw, double* u_m, double* y_m, int M, int num_mc, int num_sample)
{
    double logr, error[num_mc], error_temp, p_logr[1001];
    
    for (int ii=0; ii<1001; ii++)
    {   logr = -1.0 + ii*0.012;
        p_logr[ii] = 0;
        
        for (int j=0; j<num_mc; j++)
        {   error[j] = 0.0;
            
            for (int i=0; i<M; i++)
            {
                error_temp = y_m[i] - Flow_rate(u_m[i], rw[j], logr, Tu[j], Tl[j], L[j], Kw[j]);
                error[j] += error_temp * error_temp;
            }
            
            p_logr[ii] += likely_function(sigma_e[j], M, error[j]);
        }
    }

       Generate_sample(r_new, num_sample, p_logr, -1.0, 11.0);//Rejection method for sampling on a pdf
       cout << "New samples for r have been generated!" << endl;
    
    for(int i=0; i<num_sample; i++) r_new[i] = exp(r_new[i]);
}

void Sample_p_Tu(double* sigma_e, double* rw, double* logr, double* Tu_new, double* Tl, double* L, double* Kw, double* u_m, double* y_m, int M, int num_mc, int num_sample)
{
    double Tu, error[num_mc], error_temp, p_Tu[1001];
    
    for (int ii=0; ii<1001; ii++)
    {   Tu = (89000-3*8900) + ii*8.9*6;
        p_Tu[ii] = 0;
        
        for (int j=0; j<num_mc; j++)
        {   error[j] = 0.0;
            
            for (int i=0; i<M; i++)
            {
                error_temp = y_m[i] - Flow_rate(u_m[i], rw[j], logr[j], Tu, Tl[j], L[j], Kw[j]);
                error[j] += error_temp * error_temp;
            }
            
            p_Tu[ii] += likely_function(sigma_e[j], M, error[j]);
        }
    }

       Generate_sample(Tu_new, num_sample, p_Tu, 89000-3*8900, 89000+3*8900);//Rejection method for sampling on a pdf
       cout << "New samples for Tu have been generated!" << endl;
}

void Sample_p_Tl(double* sigma_e, double* rw, double* logr, double* Tu, double* Tl_new, double* L, double* Kw, double* u_m, double* y_m, int M, int num_mc, int num_sample)
{
    double Tl, error[num_mc], error_temp, p_Tl[1001];
    
    for (int ii=0; ii<1001; ii++)
    {   Tl = (89-3*8.9) + ii*0.0089*6;
        p_Tl[ii] = 0;
        
        for (int j=0; j<num_mc; j++)
        {   error[j] = 0.0;
            
            for (int i=0; i<M; i++)
            {
                error_temp = y_m[i] - Flow_rate(u_m[i], rw[j], logr[j], Tu[j], Tl, L[j], Kw[j]);
                error[j] += error_temp * error_temp;
            }
            
            p_Tl[ii] += likely_function(sigma_e[j], M, error[j]);
        }
    }

       Generate_sample(Tl_new, num_sample, p_Tl, 89-3*8.9, 89+3*8.9);//Rejection method for sampling on a pdf
       cout << "New samples for Tl have been generated!" << endl;
}

void Sample_p_L(double* sigma_e, double* rw, double* logr, double* Tu, double* Tl, double* L_new, double* Kw, double* u_m, double* y_m, int M, int num_mc, int num_sample)
{
    double L, error[num_mc], error_temp, p_L[1001];
    
    for (int ii=0; ii<1001; ii++)
    {   L = 1120 + ii*(1680-1120)/1000;
        p_L[ii] = 0;
        
        for (int j=0; j<num_mc; j++)
        {   error[j] = 0.0;
            
            for (int i=0; i<M; i++)
            {
                error_temp = y_m[i] - Flow_rate(u_m[i], rw[j], logr[j], Tu[j], Tl[j], L, Kw[j]);
                error[j] += error_temp * error_temp;
            }
            
            p_L[ii] += likely_function(sigma_e[j], M, error[j]);
        }
    }

       Generate_sample(L_new, num_sample, p_L, 1120.0, 1680.0);//Rejection method for sampling on a pdf
       cout << "New samples for L have been generated!" << endl;
}

void Sample_p_Kw(double* sigma_e, double* rw, double* logr, double* Tu, double* Tl, double* L, double* Kw_new, double* u_m, double* y_m, int M, int num_mc, int num_sample)
{
    double Kw, error[num_mc], error_temp, p_Kw[1001];
    
    for (int ii=0; ii<1001; ii++)
    {   Kw = 9500 + ii*(13000-9500)/1000;
        p_Kw[ii] = 0;
        
        for (int j=0; j<num_mc; j++)
        {   error[j] = 0.0;
            
            for (int i=0; i<M; i++)
            {
                error_temp = y_m[i] - Flow_rate(u_m[i], rw[j], logr[j], Tu[j], Tl[j], L[j], Kw);
                error[j] += error_temp * error_temp;
            }
            
            p_Kw[ii] += likely_function(sigma_e[j], M, error[j]);
        }
    }

       Generate_sample(Kw_new, num_sample, p_Kw, 9500.0, 13000.0);//Rejection method for sampling on a pdf
       cout << "New samples for Kw have been generated!" << endl;
}

//////////////////////////////
//////////////////////////////
int main(int argc, const char * argv[]) {

    int num_mc=10000, num_sample=10000, M=40;
    double u_m[M], y_m[M];
    Read_input(u_m, y_m);
    double rw[num_mc], logr[num_mc], Tu[num_mc], Tl[num_mc], L[num_mc], Kw[num_mc], sigma_e[num_mc];
    double rw_new[num_sample], r_new[num_sample], Tu_new[num_sample], Tl_new[num_sample], L_new[num_sample], Kw_new[num_sample], sigma_e_new[num_sample];
    Generate_random_parameter(rw, logr, Tu, Tl, L, Kw, sigma_e, num_mc);
    Sample_p_sigma_e(sigma_e_new, rw, logr, Tu, Tl, L, Kw, u_m, y_m, M, num_mc, num_sample);
    Sample_p_rw(sigma_e, rw_new, logr, Tu, Tl, L, Kw, u_m, y_m, M, num_mc, num_sample);
    Sample_p_r(sigma_e, rw, r_new, Tu, Tl, L, Kw, u_m, y_m, M, num_mc, num_sample);
    Sample_p_Tu(sigma_e, rw, logr, Tu_new, Tl, L, Kw, u_m, y_m, M, num_mc, num_sample);
    Sample_p_Tl(sigma_e, rw, logr, Tu, Tl_new, L, Kw, u_m, y_m, M, num_mc, num_sample);
    Sample_p_L(sigma_e, rw, logr, Tu, Tl, L_new, Kw, u_m, y_m, M, num_mc, num_sample);
    Sample_p_Kw(sigma_e, rw, logr, Tu, Tl, L, Kw_new, u_m, y_m, M, num_mc, num_sample);


    //write results to output file
       ofstream myfile;
       myfile.open ("Parameters based on posterior sampling.txt");
       myfile << setw(10) << left << "sigma_e"
              << setw(10) << "rw"
              << setw(10) << "r"
              << setw(10) << "Tu"
              << setw(10) << "Tl"
              << setw(10) << "L"
              << setw(10) << "Kw" << endl;
       
       for (int i=0; i<num_sample; i++)
       {
           myfile << setw(10) << left << sigma_e_new[i]
                  << setw(10) << rw_new[i]
                  << setw(10) << r_new[i]
                  << setw(10) << Tu_new[i]
                  << setw(10) << Tl_new[i]
                  << setw(10) << L_new[i]
                  << setw(10) << Kw[i] << endl;
       }
       myfile.close();
       
       cout << "Job completed! Please find results in Parameters based on posterior sampling.txt"<< endl;
       return 0;

}
