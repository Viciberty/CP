//
//  main.cpp
//  Integral
//
//  Created by RP on 5/8/17.
//  Copyright © 2017 RP. All rights reserved.
//

#include <iostream>
#include <functional>
#include <math.h>
#include "spline.h"
#include <vector>

using namespace std;
#define epsilon 0.0000001  //Bisection Accuracy
#define error 1
#define division 100000
#define hbar 6.582119513e-22


/**********************************************************
 *Interpolation
 **********************************************************/
class Interpolation
{
    int n;
    std::vector<double> x;
    std::vector<double> y;
    tk::spline s;
public:
    void _Initial(const double* x_, const double* y_, int n_)//Initialization
    {
        n = n_; //n data
        for (int i=0;i<n;i++) {x.push_back(x_[i]); y.push_back(y_[i]);}
        s.set_points(x,y);
    }
/*************************************************
*Calculate Function Value
*order<=0 3-degree spline interpolation
*order==1 linear interpolation
*order==2 quadratic interpolation
*order or more
**************************************************/
    double value(double var, int order=0)
    {
        if ((var < x[0]) || (var > x[n-1])) {cout<<"Wrong Input\n"; return 0;} //In the definition domain
        double result=0; //function value
    if (order>0) //order interpolation
    {
        int index = search(var, order, 0, n-1) / order *order;//Fine Interval
        if (index >= (n-1) - (n-1)%order) order = (n-1)%order;//The last serval points use lower orders
            double l[order+1];
            for (int i=index; i< index+order+1; i++){
                l[i-index]=1.0;
                for (int j=index; j< index+order+1; j++)
                {if (i!=j) l[i-index]*=(var - x[j])/(x[i]-x[j]);}
                result+=l[i-index]*y[i];
                }
     }
     else //3-degree spline interpolation
     {
        result = s(var);
     }
    return result;
 }


/*************************************************
*Calculate Root-Bisection
**************************************************/
//value_interval
    int* value_interval(double var, int* temp)
    {
        int count=0;
        for (int i=0; i<n; i++)
        {
            if ((y[i]-var)*(y[i+1]-var)<=0) {temp[count]=i; count++;}
            }
        return temp;
    }
//calculate root-bisection
    double root(double begin, double end, double var)
    {
        if ((end-begin)<epsilon) return (end+begin)/2;
        if ((value(begin)-var)*(value((begin+end)/2)-var) <= 0) return root(begin,(begin+end)/2,var);
         else if ((value(end)-var)*(value((begin+end)/2)-var) <= 0) return root((begin+end)/2,end,var);
         else return -1;
    }
    
    
private://search
    int search(double var, int order, int begin, int end)
    {
        for (int i=0;i<n;i++)
            if ((var >= x[i])&&(var <= x[i+1])) return i;
        return -1;
    }
};

/**********************************************************
 *Integral Function
 *order=1 trapezoid method
 *order=2 simpson method
 
 *share is the partition of the set
 **********************************************************/
double integral(double (*f_)(double,Interpolation,Interpolation,double),double begin,double end, Interpolation M, Interpolation V,double Energy, int order=2,int share = division)
{
    auto f = [&](double x){return f_(x,M,V,Energy);};
    double sum=0;
    double delta=(end-begin)*1.0/share;
    
    double a=begin;
    double b=begin,mid;
    for (int count= 0; count<share ; count++)
    {
        a = begin+count*delta;
        b = begin+(count+1)*delta;
        mid = (a+b)/2.0;
        switch (order)
        {
        case 1: {sum +=delta*f(mid); break;}
        case 2: {sum += (f(a)+ 4.0*f(mid)+f(b))*delta/6.0; break;}
        case 4: {sum += (0.0778*f(a)+0.3556*f(a+(b-a)/4)+0.1333*f(a+2*(b-a)/4)+0.3556*f(a+3*(b-a)/4)+0.0778*f(b))*delta; break;
            }
        }
    }
    return sum;
}


/**********************************************************
 *Auxiliary Function
 **********************************************************/

//Function for Assistance

double W_help(double q, Interpolation M, Interpolation V, double Energy)
{
    return sqrt(M.value(q)*(V.value(q)-Energy));
}

double T_help(double q, Interpolation M, Interpolation V, double Energy)
{
    return sqrt(M.value(q)/(Energy-V.value(q)));
}



/*******************************************************************************************************************
 *Main Program
 *******************************************************************************************************************/
int main(int argc, const char * argv[]) {
    
/**********************************************************
 *Data Initialization
 **********************************************************/
    const double q[51] = { 0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4,2.45,2.5 };
    const double V[51] = { 14.757,14.022,11.665,7.996,3.953,0.848,0,1.727,4.741,7.607,9.199,9.421,8.751,7.503,5.879,4.185,3,2.655,3.036,3.799,4.693,5.57,6.281,6.602,6.697,6.654,6.337,6.263,5.425,5.421,4.961,4.541,4.153,3.76,3.393,2.885,2.441,2.042,1.69,1.405,1.079,0.844,0.635,0.59,0.393,0.112,-0.339,-0.218,-0.329,-0.667,-1.097 };
    const double M[51] = { 271.7416,272.5246,289.6114,334.4802,424.6238,510.052,493.5916,462.7472,364.1124,246.5,275.471,327.9146,317.6718,315.2416,332.7344,323.2862,328.7846,308.0322,270.7788,272.8552,210.888,186.1278,193.7258,201.376,186.2902,189.196,157.1452,201.3586,214.687,218.3236,230.8458,228.404,227.9806,222.401,220.5856,221.241,203.5452,239.7198,261.3016,279.618,184.3646,183.4018,175.2644,173.0314,170.9782,163.2932,160.109,159.8074,157.5744,305.486,296.9339 };
    const double E_0[2] = {0.9,4.8};
/**********************************************************
*Interpolation
**********************************************************/
    Interpolation f[2];
    f[0]._Initial(q, V, 51);  //V(q)
    f[1]._Initial(q, M, 51);  //M(q)
    
 //   for (int i=0;i<251;i++) cout<<0.01*i<<" "<<f[1].value(0.01*i)<<endl;
    double root_a,root_b,root_x,root_y,root_c;
/**********************************************************
*Numerical Integration1  E_0 = 0.9
 **********************************************************/
    
    int temp_1[3];
    f[0].value_interval(0.9, temp_1);
    root_a = f[0].root(q[temp_1[0]], q[temp_1[0]+1], E_0[0]);
    root_b = f[0].root(q[temp_1[1]], q[temp_1[1]+1], E_0[0]);
    root_c = f[0].root(q[temp_1[2]], q[temp_1[2]+1], E_0[0]);
    cout<<"ROOT1: "<<root_a <<" "<<root_b<<" "<<root_c<<endl;

    double W = exp(-2*integral(W_help, root_b+epsilon, root_c-epsilon, f[1],f[0],E_0[0]));
    double T = hbar*integral(T_help, root_a+error*epsilon, root_b-error*epsilon, f[1],f[0],E_0[0]);
    
//    cout<<f[0].value(1,0)<<endl;
//    cout<<integral(W_help, root_b+error*epsilon, root_c-error*epsilon, f[1],f[0],E_0[0])<<endl;
//    cout<<integral(T_help, root_a+error*epsilon, root_b-error*epsilon, f[1],f[0],E_0[0])<<endl;
//    cout<<T_help(root_a+error*epsilon,f[1],f[0],E_0[0])*sqrt(epsilon)<<endl;
    cout<<"tau = "<<T/W<<" when E_0 = 0.9"<<endl<<endl;
 
/**********************************************************
*Numerical Integration2  E_0 = 4.8
 **********************************************************/

    int temp[5];
    f[0].value_interval(4.8, temp);
    root_a = f[0].root(q[temp[0]], q[temp[0]+1], E_0[1]);
    root_b = f[0].root(q[temp[1]], q[temp[1]+1], E_0[1]);
    root_x = f[0].root(q[temp[2]], q[temp[2]+1], E_0[1]);
    root_y = f[0].root(q[temp[3]], q[temp[3]+1], E_0[1]);
    root_c = f[0].root(q[temp[4]], q[temp[4]+1], E_0[1]);
    cout<<"ROOT2: "<< root_a <<" "<<root_b<<" "<<root_x<<" "<<root_y<<" "<<root_c<<endl;
    
    W = exp(-2*(integral(W_help, root_b+epsilon, root_x-epsilon, f[1],f[0],E_0[1])+integral(W_help, root_y+epsilon, root_c-epsilon, f[1],f[0],E_0[1])));
    T = hbar*(integral(T_help, root_a+error*epsilon, root_b-error*epsilon, f[1],f[0],E_0[1])+integral(T_help, root_x+error*epsilon, root_y-error*epsilon, f[1],f[0],E_0[1]));
//   cout<<integral(T_help, root_x, root_y, f[1],f[0],E_0[1])<<endl;
//   cout<<integral(W_help, root_b+epsilon, root_x-epsilon, f[1],f[0],E_0[1])+integral(W_help, root_y+epsilon, root_c-epsilon, f[1],f[0],E_0[1])<<" "<<(integral(T_help, root_a+error*epsilon, root_b-error*epsilon, f[1],f[0],E_0[1])+integral(T_help, root_x+error*epsilon, root_y-error*epsilon, f[1],f[0],E_0[1]))<<endl;
    
    cout<<"tau2 = "<<T/W<<" when E_0 = 4.8 "<<endl;  
 
    return 0;
}
