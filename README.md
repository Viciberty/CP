# CP
HW of Computational Physics

HW1：Numerical Integral  
####Parameter
1. 插值 double value(double var, int order=2);
  - order 
    =0 三次样条插值
    =1 线性插值
    =2 二次插值
2. 二分求根 int* value_interval(double var, int* temp);找根单调区间； double root(double begin, double end, double var);二分法查找
  - epsilon 0.000001 //二分法精度（double保证15位精度）
3. 数值积分 double integral(double (*f_)(double,Interpolation,Interpolation,double),double begin,double end, Interpolation M, Interpolation V,double Energy, int order=2,int share = division)
  - order 
    =1 梯形积分公式
    =2 Simpson积分公式  
  - error 1  //从根右边的error个epsilon开始积分
  - division 1000000//分划个数

####Environment  
System: MacOS Sierra  
Compiler: Apple LLVM version 8.1.0 (clang-802.0.42) g++  
Option: C++11  
File: main.cpp 主函数   spline.h 样条插值库  

HW2：Poisson Equation   
####Parameter  
1. 研究的区域范围  
  - r_Min r方向最小值，暂取0   
  - r_Max r方向最大值，暂取50  
  - z_Min z方向最小值，暂取-50   
  - z_Max z方向最大值，暂取50   
2. 差分方程的构建  
  - N_r r方向的划分区域个数，取的点数为格子数+1  
  - N_z z方向的划分区域个数，取的点数为格子数+1  
  - b差分方程常量部分 double b(int x)  
  - A差分方程矩阵部分 double A(int x, int y)  
3. 边界条件确定 double boundary(int i, int j, int method=1)  
  - method  
    =0 MC蒙卡确定边值条件  
    =1 多极矩展开确定边界条件  
    =2 自然边界条件（取0）    

    - MC蒙卡确定边值条件 double MC(int i,int j)  
    - 多极矩展开确定边值条件 double Quadratic(int i,int j)  
4. 迭代法确定  
  - norm确定范数 double norm(int n=-1)   
  - Jacobi方法 double Jaccobi(int M, double e=1e-3)  
    - M 最大步数  
    - e 最小误差  
  - 超松弛算法 double SOR(int M, double omega=1.407, double e=1e-3)  
    - M 最大步数  
    - e 最小误差  

####Environment  
System: MacOS Sierra  
Compiler: Apple LLVM version 8.1.0 (clang-802.0.42) g++  
Option: C++11  
File: main.cpp 主函数   
