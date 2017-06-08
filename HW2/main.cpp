//
//  main.cpp
//  PDE
//
//  Created by RP on 6/5/17.
//  Copyright © 2017 RP. All rights reserved.
//
using namespace std;
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <ctime>
#define pi 3.1415926

#define N_r 100 //r_division
#define N_z 200 //z_division
#define r_Min 0  //r_min(fm)
#define r_Max 50 //r_max(fm)
#define z_Min -50 //z_min(fm)
#define z_Max 50 //z_max(fm)
#define epsilon 1e-20 //interval_epsilon

#define step 5000 //steps
const double r_delta = r_Max*1.0/N_r;
const double z_delta = (z_Max-z_Min)*1.0/N_z;

const int r_point = N_r+1;
const int z_point = N_z+1;
double phi[r_point*z_point];
double phi_0[r_point*z_point];
double b_[r_point*z_point];
double A_[r_point*z_point][r_point*z_point];

#define z_int_Min 0
#define z_int_Max 200
#define r_int_Min 0
#define r_int_Max 100

/*****************************
 * Assistance Function
 *****************************/
double f(double r, double z)
{
    const double radius = sqrt(z*z+r*r);
    double cos = (radius==0)?0 : z / radius; //origin
    const double Y20 = 0.25*sqrt(5.0/pi)*(3*cos*cos-1);
    const double Y30 = 0.25*sqrt(7.0/pi)*(5*cos*cos*cos-3*cos);
    const double r0 = 10.0*(1+1.0*Y20+0.5*Y30);
    const double f_value = 0.8/(1+exp((radius-r0)/0.6));
    return f_value;
    }

double r(int i)
{
    if (i<r_point) {return r_Min+i*r_delta;}
    return -9999999999;
}

double z(int j)
{
    if (j<z_point){return z_Min+j*z_delta;}
    return -9999999999;
}

void print(double (*f)(double,double), double r_start=r_Min, double z_start=0, double r_n=r_point)
//print the value when z=z_start
{
    int i = 0 ;
    for (i=0;i<r_n;i++)
        cout << f(r_start+i*r_delta,z_start)<<endl;
    }


//Funtion
double MC(int i,int j);
double Quadratic(int i, int j);
double boundary(int i, int j, int method=1)
{
switch (method)
    {
    case 0:return MC(i,j);
    case 1:return Quadratic(i,j);
    case 2:return 0;
    default: return -999999999;
    }
}

double b(int x)
{
    int j = x % z_point;
    int i = x / z_point;
    if ((j==0) || (j==z_point-1) || (i==r_point-1) ) return boundary(i, j);
    return -f(r(i),z(j));
}

double A(int x, int y)
{
    int j = x % z_point;
    int i = x / z_point;
    //Inner Condition
    if ( (i>0) && (i<r_point-1) && (j>0) && (j<z_point-1))
    {
    if (y==(i-1)*z_point+j) {return 1/(r_delta*r_delta)-1/(2*r_delta*r(i));}
    if (y==(i+1)*z_point+j) {return 1/(r_delta*r_delta)+1/(2*r_delta*r(i));}
    if (y==i*z_point+j-1) {return  1/(z_delta*z_delta);}
    if (y==i*z_point+j+1) {return 1/(z_delta*z_delta);}
    if (y==i*z_point+j) {return -2*(1/(z_delta*z_delta)+1/(r_delta*r_delta));}
    }
    //r=0
    if ((i==0) && (j>0) && (j<z_point-1))
    {
        if (y==j-1) return 1.0/(z_delta*z_delta);
        if (y==j+1) return 1.0/(z_delta*z_delta);
        if (y==z_point+j) return 4.0/(r_delta*r_delta);
        if (y==j) return -(4.0/(r_delta*r_delta)+2.0/(z_delta*z_delta));
        }
    //Boundary Condition
    if ((j==0) || (j==z_point-1) || (i==r_point-1) )
        {
            if (y==i*z_point+j) return 1;
        }
    return 0;
    
}

/**************************************
 * Boundary Condition
 **************************************/
/*
void interval()//Assume an interval for integrate
{
    
    int i= 0 ;
    for (i=1;i<z_point/2;i++){
        if ((b(i)>0?b(i):-b(i))< epsilon) z_int_Min = i;
    }
    for (i=z_point/2;i<z_point-1;i++){
        if ((b(i)>0?b(i):-b(i))> epsilon) z_int_Max = i+1;
    }
    for (i=z_point/2;i<(z_point-1)*r_point;i+=z_point){
        if ((b(i)>0?b(i):-b(i))> epsilon) r_int_Max = i/z_point+1;
    }
    cout<< z_int_Min<<" "<<z_int_Max<<" "<<r_int_Max<<endl;
    cout<<z(z_int_Min)<<" "<<z(z_int_Max)<< " "<<r(r_int_Max)<<endl;
}
*/

double random(double start, double end)
{
    return start+(end-start)*rand()/(RAND_MAX );
}


double MC(int i,int j)//r(i),z(j) boudndary condition with MC
{
    const int N=1000000;
    double D = r(r_int_Max) * 2.0 * pi * (z(z_int_Max) - z(z_int_Min));
    int count;
    double sum=0;
    double xi_r,xi_z,xi_phi,distance;
    srand(unsigned(time(0)));
    for (count=0;count<N;count++)
    {
        xi_r = random(0,r(r_int_Max));
        xi_phi = random(0,2*pi);
        xi_z = random(z(z_int_Min),z(z_int_Max));
 //       cout<<xi_r<<" "<<xi_phi<<" "<<xi_z<<endl;
        distance = sqrt(xi_r*xi_r+r(i)*r(i)-2*r(i)*xi_r*cos(xi_phi)+(xi_z-z(j))*(xi_z-z(j)));
        sum+= f(xi_r,xi_z)*1.0/distance * xi_r;
        }
        return sum*1.0*D/N/4/pi;
    
    }


double Quadratic(int i,int j)//Multiple moment expansion
{
    double varphi_0,varphi_1,varphi_2;
    varphi_0 = 4551.8/sqrt(r(i)*r(i)+z(j)*z(j));
    varphi_1 = 9528.8*z(j)/(sqrt(r(i)*r(i)+z(j)*z(j))*sqrt(r(i)*r(i)+z(j)*z(j))*sqrt(r(i)*r(i)+z(j)*z(j)));
    varphi_2 = 585130*(z(j)*z(j)-0.5*r(i)*r(i))/(2*(sqrt(r(i)*r(i)+z(j)*z(j))*sqrt(r(i)*r(i)+z(j)*z(j))*sqrt(r(i)*r(i)+z(j)*z(j)))*(sqrt(r(i)*r(i)+z(j)*z(j))*sqrt(r(i)*r(i)+z(j)*z(j))*sqrt(r(i)*r(i)+z(j)*z(j)))*(sqrt(r(i)*r(i)+z(j)*z(j))*sqrt(r(i)*r(i)+z(j)*z(j))*sqrt(r(i)*r(i)+z(j)*z(j))));
    return (varphi_0+varphi_1+varphi_2)/(4*pi);
}
    
/**************************************
* Linear Equation Solution
**************************************/
double norm(int n=-1)//n-norm
{
    int count;
    double sum=0;
    for (count=0;count<z_point*r_point;count++)
    {
        sum = ((phi_0[count]-phi[count])*(phi_0[count]-phi[count])>sum)?(phi_0[count]-phi[count])*(phi_0[count]-phi[count]):sum;
    }
    return sqrt(sum);
}

double Jaccobi(int M, double e=1e-3)
{
    int count=0;
    int i,j=0;
    double sum=0;
    while (count<M)
    {
        for (i=0;i<r_point*z_point;i++){
            sum=0;
            for (j=0;j<r_point*z_point;j++)
            {
                if (i!=j) sum+=A_[i][j]*phi_0[j];
            }
            phi[i]=(b_[i]-sum)/A_[i][i];
        }
        if (norm(2)<e) {cout<<"success!"<<endl; break;}
        for (i=0;i<r_point*z_point;i++){
            phi_0[i] = phi[i];
        }
        count++;
    }
    return 0;
}

double SOR(int M, double omega=1.407, double e=1e-3)
{
    int count=0;
    int tmp,i,j=0;
    double sum=0,sum2=0;
    while (count<M)
    {
        for (i=0;i<r_point*z_point;i++){
            sum=sum2=0;
            if (i==0)
            {   for (j=1;j<r_point*z_point;j++) sum+=A_[i][j]*phi_0[j];
                phi[i]= (1-omega)*phi_0[i]+omega*(b_[i]-sum)/A_[i][i];
            }
            else if (i==r_point*z_point-1)
            {
                for (j=0;j<r_point*z_point-1;j++) sum+=A_[r_point*z_point-1][j]*phi[j];
                phi[i]= (1-omega)*phi_0[i]+omega*(b_[i]-sum)/A_[i][i];
            }
            else //其他
            {
                for (j=0;j<i;j++) sum+=A_[i][j]*phi[j];
                for (j=i+1;j<r_point*z_point;j++) sum2+=A_[i][j]*phi_0[j];
                phi[i]=(1-omega)*phi_0[i]+omega*(b_[i]-sum-sum2)/A_[i][i];
                
            }
            
        }

        if (norm(2)<e) {cout<<"success!"<<" "<<count<<endl; break;}
        for (i=0;i<r_point*z_point;i++)
        {
            phi_0[i] = phi[i];
        }
        count++;
    }
    return 0;
 }

/**************************************
 * Main Procedure
 **************************************/

int main(int argc, const char * argv[]) {
    int i,j = 0 ;

//Initialization
    for (i=0;i<r_point*z_point;i++){
        phi[i]=0;
        phi_0[i]=0;
        b_[i]=b(i);
        for (j=0;j<r_point*z_point;j++) A_[i][j]=A(i,j);
    }
 

/*
//test linear equation solution!
    for (i=0;i<r_point*z_point;i++){
        phi[i]=0;
        phi_0[i]=0;
        b_[i]=0;
        for (j=0;j<r_point*z_point;j++) A_[i][j]=0;
    }
    b_[0]=b_[1]=0;
    b_[2]=b_[3]=1;
    A_[0][0]=A_[1][1]=A_[2][2]=A_[3][3]=4;
    A_[0][1]=A_[0][2]=A_[1][0]=A_[1][3]=A_[2][0]=A_[2][3]=A_[3][1]=A_[3][2]=-1;

*/

//    interval();
//    cout<<"r="<<r(0)<<" z="<<z(200)<<" "<<MC(1,200)<<endl;
//    for (i=0;i<r_point;i++)
//        cout<<b(i*z_point+1)<<endl;
    
//Execution and Output
 //   Jaccobi(step);
    SOR(step);

    for (i=0;i<r_point;i++) 
    {
        for (j=0;j<z_point;j++) cout<<r(i)<<" "<<z(j)<<" "<<phi[i*z_point+j]<<endl;
    }

    //    cout << "Hello, World!\n";
    return 0;
}
