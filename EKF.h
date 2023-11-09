#pragma once
#include "matrix.h"
#include <ostream>
#include <vector>
#include <iostream>
#include <initializer_list>
#include "angle.h"
using namespace std;
class EKF
{
public:

	// Constructor

	EKF(int n_);                  // n: Number of estimation parameters

								  // Initialization of a new problem

	void Init(double t_, Vector x_, Matrix P_);
	void Init(double t_, Vector x_, Vector sigma);

	// Update of filter parameters

	void TimeUpdate(double        t_,    // New epoch
		const Vector& x_,    // Propagetd state 
		const Matrix& Phi);  // State transition matrix 
	void TimeUpdate(double        t_,    // New epoch
		const Vector& x_,    // Propagated state 
		const Matrix& Phi,   // State transition matrix 
		const Matrix& Qdt);  // Accumulated process noise
	void EKF::TimeUpdate(
		double        t_,
		const Vector& x_,
		const Matrix& Phi,
		const Matrix& Qdt,
		const Vector& w);
		//const Quaternion& q,
		//const Vector& b);

	void MeasUpdate(double        z,     // Measurement at new epoch
		double        g,     // Modelled measurement
		double        sigma, // Standard deviation
		const Vector& G);    // Partials dg/dx

	void MeasUpdate(const Vector& z,     // Measurement at new epoch
		const Vector& g,     // Modelled measurement
		const Vector& s,     // Standard deviation
		const Matrix& G);    // Partials dg/dx

	void MeasUpdate(const Vector& z,     // Measurement at new epoch
		const Vector& g,     // Modelled measurement
		const Matrix& Inv_W, // Measurement covariance
		const Matrix& G);    // Partials dg/dx
	void MeasUpdate(const Quaternion& z,
		const Quaternion& g,
		const Matrix& Inv_W,
		const Matrix& G);
							 // Access to filter parameters

	double   Time();                     // Time    
	Vector   State();                    // State parameters
	Matrix   Cov();                      // Covariance matrix
	Vector   StdDev();                   // Standard deviation

private:

	// Elements
	int      n;                          // Number of state parameters
	double   t;                          // Time 
	Vector   x;                          // State parameters
	Matrix   P;                          // Covariance matrix

};
// Least squares estimation class 
class LSQ
{
public:

	// Constructor
	LSQ(int nEst); // nEst: Number of estimation parameters

				   // Number of data equations
	int nData() { return n; };

	// Reset to solve a new problem
	void Init();
	void Init(const Vector& x,   // A priori parameters
		const Matrix& P);  // A priori covariance

						   // Add a data equation of form Ax=b to a least squares system
						   // (performs a row-wise QR transformation using Givens rotations)
	void Accumulate(const Vector& A, double b, double sigma = 1.0);

	// Solve the LSQ problem for vector x by backsubstitution
	void Solve(Vector& x);

	// Covariance matrix
	Matrix Cov();

	// Standard deviation
	Vector StdDev();

	// Square-root information matrix and data vector
	Matrix SRIM();      // Copy of R
	Vector Data();      // Copy of d

private:

	// Elements
	int      N;         // Number of estimation parameters
	int      n;         // Number of data equations 
	Vector   d;         // Right hand side of transformed equations
	Matrix   R;         // Square-root information matrix 
						// (Upper right triangular matrix)

};
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
void InvUpper(const Matrix& R, Matrix& T);
Vector StatePrediction(const Vector& X, const double& delta_t, const Vector& w_);
struct LowFreqErrorVector {
	vector<double> coeffs;  // 使用vector来存储系数

								 // 使用初始化列表构造函数
	LowFreqErrorVector(initializer_list<double> init = {0.0,0.0,0.0,0.0,0.0,0.0}) : coeffs(init) {
		if (coeffs.size() % 6 != 0) {
			std::cout << "输入参数维数必须为6的倍数" << std::endl;
			exit(1);  // 终止程序
		}
	}

};
inline ostream& operator<<(ostream& os, const LowFreqErrorVector& vec) {
	os << "LowFreqErrorVector: [";
	for (size_t i = 0; i < vec.coeffs.size(); ++i) {
		os << vec.coeffs[i];
		if (i < vec.coeffs.size() - 1) {
			os << ", ";
		}
	}
	os << "]";
	return os;
}
void normalize(Quaternion& q);