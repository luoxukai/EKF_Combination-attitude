#pragma once
#include <ostream>
#include "matrix.h"
#include <iostream>
#include <string>
using namespace std;


struct Quaternion   //��Ԫ��
{
	double w, x, y, z; //wΪʵ��
	//�˷�
	Quaternion operator*(double scalar) const {
		return{ w * scalar, x * scalar, y * scalar, z * scalar };
	}
	//Ĭ�Ϲ��캯��
	Quaternion() :w(1.0), x(0.0), y(0.0), z(0.0) {};
	// ���캯��
	Quaternion(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}

	// �ӷ����������
	Quaternion operator+(const Quaternion& q) const
	{
		return Quaternion(w + q.w, x + q.x, y + q.y, z + q.z);
	}

	// �������������
	Quaternion operator-(const Quaternion& q) const
	{
		return Quaternion(w - q.w, x - q.x, y - q.y, z - q.z);
	}

	// �˷����������
	Quaternion operator*(const Quaternion& q) const
	{
		double w_ = w * q.w - x * q.x - y * q.y - z * q.z;
		double x_ = w * q.x + x * q.w + y * q.z - z * q.y;
		double y_ = w * q.y - x * q.z + y * q.w + z * q.x;
		double z_ = w * q.z + x * q.y - y * q.x + z * q.w;
		return Quaternion(w_, x_, y_, z_);
	}

	// ������������أ����չ���ͷ�����
	Quaternion operator/(const Quaternion& q) const
	{
		double norm = q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
		Quaternion conj(q.w, -q.x, -q.y, -q.z);
		Quaternion result = *this * conj;
		return Quaternion(result.w / norm, result.x / norm, result.y / norm, result.z / norm);
	}
	//���
	friend ostream& operator<<(ostream& os, const Quaternion& q) {
		os << "Quaternion: [" << q.w << ", " << q.x << ", " << q.y << ", " << q.z << "]";
		return os;
	}
};
struct Vector3    //��������������
{
	double x, y, z;
	//Ĭ�Ϲ��캯��
	Vector3() : x(0.0), y(0.0), z(0.0) {};
	// ���캯��
	Vector3(double x_, double y_, double z_) :  x(x_), y(y_), z(z_) {}

};
Quaternion toA(Vector3 omega, Quaternion& omega_); //��άת��ά(x,y,z) -> (0,x,y,z)
double clamp(const double& value, const double& low, const double& high); //���������Сֵ
class cal_angle
{
public:
	static double angleBetweenQuaternions(const Quaternion& Q1, const Quaternion& Q2);//������Ԫ���н�

	static Vector3 rotate_by_quaternion(const Vector3& v, Quaternion& q); //����ͨ����Ԫ����ת

	static double angleBetweenV(const Vector3& v1, const Vector3& v2);   //���������н�

	static Quaternion conjugate(const Quaternion& V); //������Ԫ��

	static double dotProduct(const Quaternion& Q1, const Quaternion& Q2); //��Ԫ�����

	static double magnitude(const Quaternion& Q);//ģ��

	static Quaternion multiply(const Quaternion& Q1, const Quaternion& Q2); //��Ԫ������Ԫ��

	static Quaternion V2Q(const Vector3& v); //����ת��Ԫ��

	

	static double dotV(const Vector3& v1, const Vector3& v2); //�������

	static double magV(const Vector3& v); //����ģ��

	static Quaternion normalize(Quaternion& q); //��Ԫ���淶��

};
Vector V_Quaternion(Quaternion q);
struct SensorData {    //�ļ��Ľṹ��
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
// ������Ԫ������
Quaternion InverseQuaternion(const Quaternion& q);