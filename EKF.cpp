#include "EKF.h"
#include "angle.h"
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

cal_angle cal;
//------------------------------------------------------------------------------
//
// InvUpper
//
// Purpose:
//
//   Inversion of an upper right triangular matrix
//
// Input/Output:
//
//   R    Upper triangular square matrix
//   T    Inverse of R
//
// Note:
//
//   This function may be called with the same actual parameter for R and T
//
//------------------------------------------------------------------------------

void InvUpper(const Matrix& R, Matrix& T)
{

	const int N = R.size1();   // Dimension of R and T

	int    i, j, k;
	double Sum;

	if (R.size2() != N || T.size1() != N || T.size2() != N) {
		cerr << " ERROR: Incompatible shapes in InvUpper" << endl;
		return;
	}

	// Check diagonal elements

	for (i = 0; i<N; i++)
		if (R(i, i) == 0.0) {
			cerr << " ERROR: Singular matrix in InvUpper" << endl;
			return;
		}
		else {
			// Compute the inverse of i-th diagonal element.
			T(i, i) = 1.0 / R(i, i);
		};

		// Calculate the inverse T = R^(-1)

		for (i = 0; i<N - 1; i++)
			for (j = i + 1; j<N; j++) {
				Sum = 0.0; for (k = i; k <= j - 1; k++)  Sum += T(i, k)*R(k, j);
				T(i, j) = -T(j, j)*Sum;
			};


}

//------------------------------------------------------------------------------
//
// StatePrediction
//
// Purpose:
//
//   计算状态量的一步预测
//
// Input/Output:
//
//   X    状态量
//   X_   状态量一步预测
//
//------------------------------------------------------------------------------




Vector StatePrediction(const Vector& X, const double& delta_t,const Vector& w)
{
	//姿态误差
	//cout << "X:" << X << endl;
	Vector delta_q(X(0), X(1),X(2));
	//q.w = X(0), q.x = X(1), q.y = X(2), q.z = X(3);
	//陀螺仪常漂误差
	//Vector delta_b(X(3), X(4), X(5));
	//四维表示
	//Quaternion A_delta_b,A_omega_;
	//toA(delta_b, A_delta_b);
	//toA(omega_, A_omega_);
	//低频噪声系数误差
	//LowFreqErrorVector k{ X(6), X(7), X(8), X(9),X(10),X(11) };
	//double k_ax1 = X(7), k_bx1 = X(8);
	//double k_ay1 = X(9), k_by1 = X(10);
	//double k_az1 = X(11), k_bz1 = X(12);
	//陀螺角速度估计值
	//Vector w_ = w - delta_b;
	Vector w_ = w;
	//姿态误差估计
	//Quaternion delta_q_ = delta_q+((delta_q*A_delta_b) *0.5 + (delta_q*A_omega_)*0.5 - (A_omega_*delta_q)*0.5)*delta_t;
	//角速度估计值反对称矩阵
	Matrix wx(3,3);
	wx.M[0][1] = w_(2), wx.M[0][2] = -w_(1),
		wx.M[1][0] = -w_(2), wx.M[1][2] = w_(0),
		wx.M[2][0] = w_(1), wx.M[2][1] = -w_(0);
	//姿态误差变化率
	Vector delta_q_rate = wx*delta_q;
	//姿态误差一步预测
	Vector pre_delta_q = delta_q + delta_q_rate*delta_t;
	//cout << "delta_q.w:"<<delta_q_.w << endl;
	//陀螺常漂误差不变
	//低频噪声系数误差不变
	/*
	//姿态估计
	Quaternion pre_q;
	pre_q = q + A_omega_*q*0.5*delta_t;
	//陀螺常漂估计
	Vector pre_b;
	pre_b = b;
	//低频噪声系数估计
	LowFreqErrorVector pre_k;
	pre_k = k;
	*/
	//Vector result;
	Vector X_(X.size());
	X_(0) = pre_delta_q(0), X_(1) = pre_delta_q(1), X_(2) = pre_delta_q(2);
	//X_(3) = delta_b(0), X_(4) = delta_b(1), X_(5) = delta_b(2);
	//X_(6) = k.coeffs[0], X_(7) = k.coeffs[1], X_(8) = k.coeffs[2], X_(9) = k.coeffs[3], X_(10) = k.coeffs[4], X_(11) = k.coeffs[5];
	//cout << "X_:" << X_ << endl;
	return X_;
}





									   //------------------------------------------------------------------------------
									   //
									   // EKF class (implementation)
									   //
									   //------------------------------------------------------------------------------

									   //
									   // Constructor
									   //

EKF::EKF(int n_)
	: n(n_),             // Number of estimation parameters
	t(0.0)             // Epoch
{
	// Allocate storage and initialize to zero
	x = Vector(n);
	P = Matrix(n, n);
}

//
// Initialization of a new problem
//

void EKF::Init(double t_, Vector x_, Matrix P_)
{
	t = t_; x = x_; P = P_;
}

void EKF::Init(double t_, Vector x_, Vector sigma)
{
	t = t_; x = x_;
	P = 0.0;
	for (int i = 0; i<n; i++) P(i, i) = sigma(i)*sigma(i);
}


//
// Access to filter parameters
//

double EKF::Time() { return t; };     // Time

Vector EKF::State() { return x; };     // State parameters

Matrix EKF::Cov() { return P; };     // Covariance matrix

Vector EKF::StdDev() {                 // Standard deviation
	Vector Sigma(n);
	for (int i = 0; i<n; i++) Sigma(i) = sqrt(P(i, i));
	return Sigma;
}

//
// Time Update
//

void EKF::TimeUpdate(
	double        t_,
	const Vector& x_,
	const Matrix& Phi,
	const Matrix& Qdt,
	const Vector& w)
	//const Quaternion& q,
	//const Vector& b)
{
	//cout << "Phi:" << Phi << endl;
	//cout << "x_:" << x_ << endl;
	t = t_;                          // Next time step
									 //x = x_;                          // Propagated state
	x = StatePrediction(x_, t_, w);
	P = Phi*P*Transp(Phi) + Qdt;     // Propagated covariance + noise
}

//
// Scalar Measurement Update
//

//
// Vector Measurement Update
//

void EKF::MeasUpdate(const Vector& z,
	const Vector& g,
	const Matrix& Inv_W,
	const Matrix& G)
{

	Matrix K(n, z.size());                 // Kalman gain

										   // Kalman gain

	K = P*Transp(G)*Inv(Inv_W + G*P*Transp(G));
	//cout << "K:" << K << endl;
	// State update
	Quaternion Z1(z(0), z(1), z(2), z(3));
	Quaternion G1(g(0), g(1), g(2), g(3));
	Quaternion Z1_G1;
	Z1_G1 = Z1 / G1;
	Vector z_g(4);
	
	z_g(0) = Z1_G1.w, z_g(1) = Z1_G1.x, z_g(2) = Z1_G1.y, z_g(3) = Z1_G1.z;
	//cout << "新息：" << Z1_G1 << endl;
	x = x + K*(z_g);
	//cout << "乘K之后的状态量：" << x << endl;
	// Covariance update

	P = (Id(n) - K*G)*P;

}

void EKF::MeasUpdate(const Quaternion& z,
	const Quaternion& g,
	const Matrix& Inv_W,
	const Matrix& G)
{

	Matrix K(n, 3);                 // Kalman gain

										   // Kalman gain

	K = P*Transp(G)*Inv(Inv_W + G*P*Transp(G));
	//cout << "K:" << K << endl;
	// State update
	Quaternion Z1_G1;
	Z1_G1 = InverseQuaternion(g)*z;
	normalize(Z1_G1);
	Vector z_g(3);
	//cout << "乘K之前的状态量q：" << Z1_G1<< endl;
	z_g(0) = Z1_G1.x, z_g(1) = Z1_G1.y, z_g(2) = Z1_G1.z;
	//cout << "新息：" << Z1_G1 << endl;
	
	x = K*(z_g);
	//Quaternion temp(1,x[1],x[2],x[3]);
	//normalize(temp);
	//x[0] = temp.w, x[1] = temp.x, x[2] = temp.y, x[3] = temp.z;
	//cout << "乘K之后的状态量q：" << x[0]<<","<<x[1] << "," << x[2] <<  endl;
	//cout << "乘K之后的状态量b：" << x[3] << "," << x[4] << "," << x[5]  << endl;
	//cout << "乘K之后的状态量k：" << x[6] << "," << x[7] << "," << x[8] << "," << x[9] << "," << x[10] << "," << x[11] << endl;
	// Covariance update

	P = (Id(n) - K*G)*P*Transp(Id(n) - K*G)+K*Inv_W*Transp(K);
	//cout << "p" << P << endl;
}
void normalize(Quaternion& q) {
	//std::cout << "Before normalization: " << q << std::endl;

	double norm = sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
	if (norm > 0) { // 避免除以零
		q.w /= norm;
		q.x /= norm;
		q.y /= norm;
		q.z /= norm;
	}
	else {
		// 无法规范化零四元数，可能需要设置为默认单位四元数或抛出错误
	}
	//std::cout << "After normalization: " << q << std::endl;
}