#include "matrix.h"
#define Mat math::matrix
template<int M, int N>  //M:输入维度，N：输出维度
class rls   
{
private:
Mat<M,1> k = {0.0}; //rls_gain
Mat<M, M> I = Mat<M, M>::identity();
float alpha = 1.0;  
Mat<M,M>P = alpha * I; //协方差类矩阵
Mat<M,N> theta = {0.0}; //参数向量

Mat<N,1> error = 0.0; //预测误差
float lambda = 0.99; //遗忘因子,默认为0.99

public:
   Mat<N,1> y = {0.0}; //当前测量值
   Mat<M,1> x = {0.0}; //当前输入
   rls(Mat<N,1>y_,Mat<M,1>x_):y(y_),x(x_){};
   rls():y(Mat<N,1>::zeros()),x(Mat<M,1>::zeros()){};
   void update(){
   const Mat<M,N>lst_theta = theta;
   const Mat<M,M>lst_P = P;
   const Mat<N,1>y_hat = lst_theta.T * x;

      error = y - y_hat;
      k = (lst_P * x)/ (lambda + x.T * lst_P * x);
      theta = lst_theta + k * error.T;
      P = (1.0f/lambda) * (lst_P - k * x.T * lst_P);
   }
};


