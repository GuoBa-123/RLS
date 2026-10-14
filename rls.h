#pragma once

#include "matrix.h"
#define MAX_P_TRACE 1e6f  //最大
#define Mat math::matrix
template<int M, int N>  //M:输入维度，N：输出维度
class rls
{
private:
   Mat<M,1> k{};          //rls_gain
   float alpha = 1000.0f; //P 初值系数 P0 = alpha * I
   Mat<M,M> P = Mat<M,M>::identity() * alpha; //协方差类矩阵
   Mat<N,1> error{};      //预测误差（先验）
   float lambda = 0.99f;  //遗忘因子,默认为0.99

public:
   typedef struct
   {
      Mat<N,1> y{};      //当前测量值
      Mat<M,1> x{};      //当前输入
      Mat<M,N> theta{};  //参数向量
   }rls_t;
   rls_t state;

   rls() = default;
   rls(float lambda_) : lambda(lambda_) { reset(); }     //只输入遗忘因子，alpha默认1000
   rls(float lambda_, float alpha_) : alpha(alpha_), lambda(lambda_) { reset(); }

   void reset() {
      P = Mat<M,M>::identity() * alpha;
      state.theta = Mat<M,N>::zeros();
      k = Mat<M,1>::zeros();
      error = Mat<N,1>::zeros();
   }
   void setLambda(float l) { lambda = l; }

   void update(const Mat<M,1> &input, const Mat<N,1> &output) {
      state.x = input;
      state.y = output;
      calculate();
   }
   void calculate(){
      const Mat<M,1> Px = P * state.x;  // O(M²)
      const float denom = lambda + (state.x.T() * Px)(0, 0);
      if (denom < 1e-6f) return;        // 数值保护

      error = state.y - state.theta.T() * state.x;    //y - y_hat
      k = Px / denom;
      state.theta = state.theta + k * error.T();
      P = (P - k * Px.T()) * (1.0f / lambda); // 利用 P 对称，免去 M³
      P = (P + P.T()) * 0.5f;                 // 保持对称
      const float tr = P.trace();
      if (tr > MAX_P_TRACE) P = P * (MAX_P_TRACE / tr); // 防 windup
   }
   const rls_t* getRlsMessage()const{
      return &state;
   }
   
};
