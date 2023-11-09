#pragma once
#include <ostream>
#include "matrix.h"
#include <iostream>
#include <string>
using namespace std;


struct Quaternion   //四元数
{
	double w, x, y, z; //w为实部
	//乘法
	Quaternion operator*(double scalar) const {
		return{ w * scalar, x * scalar, y * scalar, z * scalar };
	}
	//默认构造函数
	Quaternion() :w(1.0), x(0.0), y(0.0), z(0.0) {};
	// 构造函数
	Quaternion(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}

	// 加法运算符重载
	Quaternion operator+(const Quaternion& q) const
	{
		return Quaternion(w + q.w, x + q.x, y + q.y, z + q.z);
	}

	// 减法运算符重载
	Quaternion operator-(const Quaternion& q) const
	{
		return Quaternion(w - q.w, x - q.x, y - q.y, z - q.z);
	}

	// 乘法运算符重载
	Quaternion operator*(const Quaternion& q) const
	{
		double w_ = w * q.w - x * q.x - y * q.y - z * q.z;
		double x_ = w * q.x + x * q.w + y * q.z - z * q.y;
		double y_ = w * q.y - x * q.z + y * q.w + z * q.x;
		double z_ = w * q.z + x * q.y - y * q.x + z * q.w;
		return Quaternion(w_, x_, y_, z_);
	}

	// 除法运算符重载（按照共轭和范数）
	Quaternion operator/(const Quaternion& q) const
	{
		double norm = q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
		Quaternion conj(q.w, -q.x, -q.y, -q.z);
		Quaternion result = *this * conj;
		return Quaternion(result.w / norm, result.x / norm, result.y / norm, result.z / norm);
	}
	//输出
	friend ostream& operator<<(ostream& os, const Quaternion& q) {
		os << "Quaternion: [" << q.w << ", " << q.x << ", " << q.y << ", " << q.z << "]";
		return os;
	}
};
struct Vector3    //向量的三个分量
{
	double x, y, z;
	//默认构造函数
	Vector3() : x(0.0), y(0.0), z(0.0) {};
	// 构造函数
	Vector3(double x_, double y_, double z_) :  x(x_), y(y_), z(z_) {}

};
Quaternion toA(Vector3 omega, Quaternion& omega_); //三维转四维(x,y,z) -> (0,x,y,z)
double clamp(const double& value, const double& low, const double& high); //限制最大最小值
class cal_angle
{
public:
	static double angleBetweenQuaternions(const Quaternion& Q1, const Quaternion& Q2);//计算四元数夹角

	static Vector3 rotate_by_quaternion(const Vector3& v, Quaternion& q); //向量通过四元数旋转

	static double angleBetweenV(const Vector3& v1, const Vector3& v2);   //计算向量夹角

	static Quaternion conjugate(const Quaternion& V); //共轭四元数

	static double dotProduct(const Quaternion& Q1, const Quaternion& Q2); //四元数点积

	static double magnitude(const Quaternion& Q);//模长

	static Quaternion multiply(const Quaternion& Q1, const Quaternion& Q2); //四元数乘四元数

	static Quaternion V2Q(const Vector3& v); //向量转四元数

	

	static double dotV(const Vector3& v1, const Vector3& v2); //向量点积

	static double magV(const Vector3& v); //向量模长

	static Quaternion normalize(Quaternion& q); //四元数规范化

};
Vector V_Quaternion(Quaternion q);
struct SensorData {    //文件的结构体
	string timestamp;
	//double StarA_q0, StarA_q1, StarA_q2, StarA_q3;
	//double StarB_q0, StarB_q1, StarB_q2, StarB_q3;
	//double GyroA_X, GyroA_Y, GyroA_Z, GyroA_T;
	//double GyroB_X, GyroB_Y, GyroB_Z, GyroB_T;
	double StarA_q0 = 0.0, StarA_q1 = 0.0, StarA_q2 = 0.0, StarA_q3 = 0.0;
	double StarB_q0 = 0.0, StarB_q1 = 0.0, StarB_q2 = 0.0, StarB_q3 = 0.0;
	double GyroA_X = 0.0, GyroA_Y = 0.0, GyroA_Z = 0.0, GyroA_T = 0.0;
	double GyroB_X = 0.0, GyroB_Y = 0.0, GyroB_Z = 0.0, GyroB_T = 0.0;

	friend ostream& operator<<(std::ostream& os, const SensorData& data) {
		os << "Timestamp: " << data.timestamp << std::endl;
		os << "StarA Quaternion: [" << data.StarA_q0 << ", " << data.StarA_q1 << ", " << data.StarA_q2 << ", " << data.StarA_q3 << "]" << std::endl;
		os << "StarB Quaternion: [" << data.StarB_q0 << ", " << data.StarB_q1 << ", " << data.StarB_q2 << ", " << data.StarB_q3 << "]" << std::endl;
		os << "GyroA: [" << data.GyroA_X << ", " << data.GyroA_Y << ", " << data.GyroA_Z << ", " << data.GyroA_T << "]" << std::endl;
		os << "GyroB: [" << data.GyroB_X << ", " << data.GyroB_Y << ", " << data.GyroB_Z << ", " << data.GyroB_T << "]";
		return os;
	}
};
struct Euler {
	double roll, pitch, yaw;
};

Euler quaternionToEuler(Quaternion q);
// 计算四元数的逆
Quaternion InverseQuaternion(const Quaternion& q);