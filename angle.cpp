#include <cmath>
#include "angle.h"
#include <iostream>
#include <algorithm>
#include "matrix.h"
Vector V_Quaternion(Quaternion q) 
{
	Vector Q;
	Q(0) = q.x, Q(1) = q.y, Q(0) = q.z;
	return Q;
}
Quaternion toA(Vector3 omega, Quaternion& omega_) //��άת��ά(x,y,z) -> (0,x,y,z)
{
	omega_.w = 0;
	omega_.x = omega.x;
	omega_.y = omega.y;
	omega_.z = omega.z;
	return omega_;
}
double cal_angle::dotProduct(const Quaternion& Q1, const Quaternion& Q2) //��Ԫ�����
{
	return Q1.w*Q2.w + Q1.x*Q2.x + Q1.y*Q2.y + Q1.z*Q2.z;
}
double cal_angle::magnitude(const Quaternion& Q) //ģ��
{
	return std::sqrt(Q.w*Q.w + Q.x*Q.x + Q.y*Q.y + Q.z*Q.z);
}
double cal_angle::angleBetweenQuaternions(const Quaternion& Q1, const Quaternion& Q2)//������Ԫ���н�
{
	double dot = dotProduct(Q1, Q2);
	double mag1 = magnitude(Q1);
	double mag2 = magnitude(Q2);
	double angle = std::acos(dot / (mag1*mag2));
	return angle*(180.0/3.1415926535);
}

Quaternion cal_angle::multiply(const Quaternion& Q1, const Quaternion& Q2) //��Ԫ������Ԫ��
{
	Quaternion ans;
	double d1, d2, d3, d4;

	d1 = Q1.w * Q2.w;
	d2 = -Q1.x * Q2.x;
	d3 = -Q1.y * Q2.y;
	d4 = -Q1.z * Q2.z;
	ans.w = d1 + d2 + d3 + d4;

	d1 = Q1.w * Q2.x;
	d2 = Q2.w * Q1.x;
	d3 = Q1.y * Q2.z;
	d4 = -Q1.z * Q2.y;
	ans.x = d1 + d2 + d3 + d4;

	d1 = Q1.w * Q2.y;
	d2 = Q2.w * Q1.y;
	d3 = Q1.z * Q2.x;
	d4 = -Q1.x * Q2.z;
	ans.y = d1 + d2 + d3 + d4;

	d1 = Q1.w * Q2.z;
	d2 = Q2.w * Q1.z;
	d3 = Q1.x * Q2.y;
	d4 = -Q1.y * Q2.x;
	ans.z = d1 + d2 + d3 + d4;

	return ans;
}
Quaternion cal_angle::V2Q(const Vector3& v) //����ת��Ԫ��
{
	Quaternion q;
	q.w = 0.0;
	q.x = v.x;
	q.y = v.y;
	q.z = v.z;
	return q;
}

Quaternion cal_angle::conjugate(const Quaternion& V)//������Ԫ��
{
	Quaternion ans;
	ans.w = V.w;
	ans.x = -V.x;
	ans.y = -V.y;
	ans.z = -V.z;
	return ans;
}

Vector3 cal_angle::rotate_by_quaternion(const Vector3& v, Quaternion& q)
{
	Quaternion v_qua = V2Q(v);                                     //����תΪ������Ԫ��
	Quaternion q1 = normalize(q);   //��λ��
	Quaternion q_conjugate = conjugate(q1);                         //���㹲��
	Quaternion result = multiply(multiply(q1, v_qua), q_conjugate); //ִ����ת
	return{ result.x, result.y, result.z };                        //��ȡ��ת����鲿
}

double cal_angle::dotV(const Vector3& v1, const Vector3& v2) //�������
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}
double cal_angle::magV(const Vector3& v) //����ģ��
{
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

double clamp(const double& value, const double& low, const double& high) {
	return std::min(std::max(value, low), high);
}

double cal_angle::angleBetweenV(const Vector3& v1, const Vector3& v2)
{
	double dot = dotV(v1, v2);
	double mag1 = magV(v1);
	double mag2 = magV(v2);
	//std::cout << dot / (mag1 * mag2) << std::endl;
	double cos_theta =clamp( dot / (mag1 * mag2),-1.0,1.0);       //����cos������Сֵ
	return acos(cos_theta)*(180.0 / 3.1415926535);
}

Quaternion cal_angle::normalize(Quaternion& q)
{
	double mag =  magnitude(q);
	if (mag == 0.0) {
		std::cout << "������Ԫ��ģ��Ϊ0�����!!!!!" << std::endl;
	}
	Quaternion normalizeQ;
	normalizeQ.w = q.w / mag;
	normalizeQ.x = q.x / mag;
	normalizeQ.y = q.y / mag;
	normalizeQ.z = q.z / mag;

	return normalizeQ;
}
//��Ԫ��ת��ΪXYZ˳���ŷ����
Euler quaternionToEuler(Quaternion q) {
	Euler e;

	// Roll (X-axis rotation)
	double sinr_cosp = 2 * (q.w * q.z + q.y * q.x);
	double cosr_cosp = 1 - 2 * (q.z * q.z + q.y * q.y);
	e.roll = std::atan2(sinr_cosp, cosr_cosp)*(180.0 / 3.14159265358979323846);

	// Pitch (Y-axis rotation)
	double sinp = 2 * (q.w * q.y - q.z * q.x);
	if (std::abs(sinp) >= 1) {
		e.pitch = std::copysign(3.14159265358979323846 / 2, sinp)*(180.0 / 3.14159265358979323846); // Use 90 degrees if out of range
	}
	else {
		e.pitch = std::asin(sinp)*(180.0 / 3.14159265358979323846);
	}

	// Yaw (Z-axis rotation)
	double siny_cosp = 2 * (q.w * q.x + q.z * q.y);
	double cosy_cosp = 1 - 2 * (q.y * q.y + q.x * q.x);
	e.yaw = std::atan2(siny_cosp, cosy_cosp)*(180.0 / 3.14159265358979323846);

	return e;
}
// ������Ԫ������
Quaternion InverseQuaternion(const Quaternion& q) {
	double norm = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
	Quaternion inverse;

	// �����Ԫ����ģ���ӽ��㣬�޷�������
	if (norm < 1e-8) {
		std::cerr << "Error: Quaternion has near-zero norm, inverse does not exist." << std::endl;
		// ����һ����Ч����Ԫ��
		return{ 0, 0, 0, 0 };
	}

	double invNorm = 1.0 / norm;
	inverse.w = q.w * invNorm;
	inverse.x = -q.x * invNorm;
	inverse.y = -q.y * invNorm;
	inverse.z = -q.z * invNorm;

	return inverse;
}

