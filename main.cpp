#include "angle.h"
#include<cstdlib>
#include<iostream>
#include<fstream>
//#include<ostream>
#include<string>
#include <vector>
#include <sstream>
#include <iomanip> // 包含setprecision和fixed
#include <chrono>
#include "EKF.h"
#include "matrix.h"
using namespace std;
//解析时间戳
tm parse_timestamp(const string& timestamp) {
	tm tm;
	istringstream ss(timestamp.substr(0,19));
	ss >> get_time(&tm, "%Y-%m-%d %H:%M:%S");
	return tm;
}
int parse_milliseconds(const string& timestamp) {
	return stoi(timestamp.substr(20, 4));
}

#if 1
int main()
{


	ifstream file("test_data.txt");   //读取文件
	ofstream file_out; //输出文件
	int previous_day = -1; //初值，保证每天都有输出
	if (!file)                  //判断能否打开
	{
		cerr << "无法打开文件！" << endl;
		return 1;
	}
	string line;
	string header_line;
	getline(file, header_line); //忽略标题行
	vector<SensorData> data;
	SensorData entry;
	while (getline(file, line))//读取数据
	{
		istringstream iss(line);


		string date, time;
		iss >> date >> time;
		entry.timestamp = date + " " + time;
		iss >> entry.StarA_q0 >> entry.StarA_q1 >> entry.StarA_q2 >> entry.StarA_q3;
		iss >> entry.StarB_q0 >> entry.StarB_q1 >> entry.StarB_q2 >> entry.StarB_q3;
		iss >> entry.GyroA_X >> entry.GyroA_Y >> entry.GyroA_Z >> entry.GyroA_T;
		iss >> entry.GyroB_X >> entry.GyroB_Y >> entry.GyroB_Z >> entry.GyroB_T;


		if (entry.StarA_q0 != 0 && entry.StarA_q1 != 0 && entry.StarA_q2 != 0 && entry.StarA_q3 != 0 &&
			entry.StarB_q0 != 0 && entry.StarB_q1 != 0 && entry.StarB_q2 != 0 && entry.StarB_q3 != 0 && entry.GyroA_Y != 0)
		{
			data.push_back(entry);
		}

	}

	file.close();

	for (int i = 0; i < data.size() - 1; ++i)
	{
		SensorData& currentEntry = data[i];
		SensorData& nextEntry = data[i + 1];
		Vector currentQB(4);
		Vector nextQB(4);
		currentQB(0) = currentEntry.StarB_q0, currentQB(1) = currentEntry.StarB_q1, currentQB(2) = currentEntry.StarB_q2, currentQB(3) = currentEntry.StarB_q3;
		nextQB(0) = nextEntry.StarB_q0, nextQB(1) = nextEntry.StarB_q1, nextQB(2) = nextEntry.StarB_q2, nextQB(3) = nextEntry.StarB_q3;
		Vector diffb = nextQB - currentQB;
		Vector currentQA(4);
		Vector nextQA(4);
		currentQA(0) = currentEntry.StarA_q0, currentQA(1) = currentEntry.StarA_q1, currentQA(2) = currentEntry.StarA_q2, currentQA(3) = currentEntry.StarA_q3;
		nextQA(0) = nextEntry.StarA_q0, nextQA(1) = nextEntry.StarA_q1, nextQA(2) = nextEntry.StarA_q2, nextQA(3) = nextEntry.StarA_q3;
		Vector diffa = nextQA - currentQA;
		double sqB = sqrt(diffb(0)*diffb(0) + diffb(1)*diffb(1) + diffb(2)*diffb(2) + diffb(3)*diffb(3));
		double sqA = sqrt(diffa(0)*diffa(0) + diffa(1)*diffa(1) + diffa(2)*diffa(2) + diffa(3)*diffa(3));
		if (sqA > 1.5)
		{
			nextEntry.StarA_q0 = -nextEntry.StarA_q0;
			nextEntry.StarA_q1 = -nextEntry.StarA_q1;
			nextEntry.StarA_q2 = -nextEntry.StarA_q2;
			nextEntry.StarA_q3 = -nextEntry.StarA_q3;
		}
		if (sqB > 1.5)
		{
			nextEntry.StarB_q0 = -nextEntry.StarB_q0;
			nextEntry.StarB_q1 = -nextEntry.StarB_q1;
			nextEntry.StarB_q2 = -nextEntry.StarB_q2;
			nextEntry.StarB_q3 = -nextEntry.StarB_q3;
		}

	}
	int n = 3;
	EKF ekf(n);
	Vector tr_x0(n);//初始化状态量
	Matrix P0(1e-3, n, n); //初始协方差矩阵
	//cout << P0 << endl;
	//double t = 0; //初始时间
	Matrix Phi(n, n); //状态转移矩阵
	Matrix I(1, n, n);
	//cout << "Phi:" << Phi << endl;
	Matrix Q(1e-6, n, n); //过程噪声协方差
	//Q.M[3][3] = 1e-10, Q.M[4][4] = 1e-10, Q.M[5][5] = 1e-10;
	//Q.M[6][6] = 0.002, Q.M[7][7] = 0.002, Q.M[8][8] = 0.002, Q.M[9][9] = 0.002, Q.M[10][10] = 0.002, Q.M[11][11] = 0.002;
	Matrix R(1e-12, 3, 3); //量测噪声协方差
	//R.M[1][1] = 1e-6;
	Matrix F(n, n); //变化率状态转移矩阵
	//Vector b(3);//陀螺常漂
	//b(0) = 0, b(1) = 0, b(2) = 0;
	//Vector x_(n);   //一步转移估计值（姿态四元数，陀螺常漂，低频噪声系数）
	tm prev_tm = {};//初始时刻的时间（不包括毫秒）
	//Quaternion z(0.2368387363, 0.6148434139, 0.2842621221, -0.6964696863); //姿态初值B
	Quaternion z(0.6686909685, 0.4073784283, -0.0343741168, 0.6210584538); //A
	normalize(z);
	//std::cout << "aa:" << z << endl;
	prev_tm.tm_year = 2023 - 1900;  // 年份从1900开始
	prev_tm.tm_mon = 7 - 1;  // 月份从0开始
	prev_tm.tm_mday = 8;  // 日期
	prev_tm.tm_hour = 0;  // 小时
	prev_tm.tm_min = 0;  // 分钟
	prev_tm.tm_sec = 0;  // 秒
	int prev_ms = 105;//初始时刻的毫秒
	vector<double> timeIntervals;
#if 0	//保存时间间隔
	ofstream outFile("output_dt2.dat");
	if (!outFile) {
		std::cerr << "无法打开或创建输出文件！" << std::endl;
		return 1;
	}
	for (const auto& entry : data)
	{
		//cout << entry.timestamp << endl;
		//outFile << entry.timestamp;  // 使用重载的 operator<<
		//outFile << std::endl;
		//cout << entry.timestamp.substr(20, 3) << endl;
		tm curr_tm = parse_timestamp(entry.timestamp.substr(0, 19));
		int curr_ms = parse_milliseconds(entry.timestamp);
		time_t prev_time = mktime(&prev_tm);
		time_t curr_time = mktime(&curr_tm);
		double dt = difftime(curr_time, prev_time) + (curr_ms - prev_ms) / 1000.0;
		outFile << dt;
		outFile << endl;
		prev_tm = curr_tm;
		prev_ms = curr_ms;
	}
#endif

	
	//Quaternion a(0, 0, 1e-5, 0);
	//Quaternion b = z*a;
	//储存每一刻的估计值和协方差阵
	//vector<Vector> stateVectorHistory_f;
	//vector<Matrix> covarianceMatrixHistory_f;
	//vector<Vector> stateVectorHistory_b;
	//vector<Matrix> covarianceMatrixHistory_b;
	double W = 0.001; //卫星运动角速度，弧度/秒
	ekf.Init(0.0, tr_x0, P0); //使用初始状态和协方差初始化EKF
	Euler e_after;//欧拉角表示误差
	Euler e_before;
	Euler e_q_;//姿态估计值欧拉角
	//星敏感器安装矩阵
	//B
	//Quaternion q_install(0.214078959247460, 0.455497380651152, 0.783345592368603, -0.364778109904485);
	//C-A
	Quaternion q_install(0.211415443800940, 0.450275982216840, 0.786500835043359, -0.366019973186142);
	//C-B
	//Quaternion q_install(-0.110071454696354, 0.874480630010094, -0.238739105588471, 0.407642016819377);
	Matrix Ontology2inertial_q(3, 3);
	Ontology2inertial_q.M[0][0] = 1 - 2 * (q_install.y*q_install.y + q_install.z*q_install.z), Ontology2inertial_q.M[0][1] = 2 * (q_install.x*q_install.y - q_install.z*q_install.w), Ontology2inertial_q.M[0][2] = 2 * (q_install.x*q_install.z + q_install.y*q_install.w);
	Ontology2inertial_q.M[1][0] = 2 * (q_install.x*q_install.y + q_install.w*q_install.z), Ontology2inertial_q.M[1][1] = 1 - 2 * (q_install.x*q_install.x + q_install.z*q_install.z), Ontology2inertial_q.M[1][2] = 2 * (q_install.y*q_install.z - q_install.x*q_install.w);
	Ontology2inertial_q.M[2][0] = 2 * (q_install.x*q_install.z - q_install.y*q_install.w), Ontology2inertial_q.M[2][1] = 2 * (q_install.y *q_install.z + q_install.x*q_install.w), Ontology2inertial_q.M[2][2] = 1 - 2 * (q_install.x*q_install.x + q_install.y*q_install.y);

#if 1
	for (const auto& entry : data) {
		//cout << "P:" << P0 << endl;
		//cout << "姿态初值：" << z;
		//cout << entry << endl;
		//tm tm = parse_timestamp(entry.timestamp.substr(0, 23));
		tm curr_tm = parse_timestamp(entry.timestamp.substr(0, 19));
		int curr_ms = parse_milliseconds(entry.timestamp);
		time_t prev_time = mktime(&prev_tm);
		time_t curr_time = mktime(&curr_tm);
		//cout << entry.timestamp << endl;
		double dt = difftime(curr_time, prev_time) + (curr_ms - prev_ms) / 1000.0;  //当前时刻和前一时刻的时间间隔
		//量测模型估计值
		Vector g_(3);
		//星敏量测值输入
		//Vector g_input(4);
		//g_input(0) = entry.StarA_q0, g_input(1) = entry.StarA_q1, g_input(2) = entry.StarA_q2, g_input(3) = entry.StarA_q3;
		//	g_input = Ontology2inertial_q*g_input;
		//量测值g
		Quaternion g(entry.StarA_q0, entry.StarA_q1, entry.StarA_q2, entry.StarA_q3);
		g = g*q_install;
		normalize(g);
		//陀螺载体系到本体系旋转矩阵
		Matrix gyr2ntology(3, 3);
		gyr2ntology.M[0][0] = -1, gyr2ntology.M[0][1] = 0, gyr2ntology.M[0][2] = 0;
		gyr2ntology.M[1][0] = 0, gyr2ntology.M[1][1] = -1, gyr2ntology.M[1][2] = 0;
		gyr2ntology.M[2][0] = 0, gyr2ntology.M[2][1] = 0, gyr2ntology.M[2][2] = 1;

		Vector w_input(entry.GyroA_X, entry.GyroA_Y, entry.GyroA_Z);

		w_input = gyr2ntology*(w_input*(3.141592653589793/180.0));
		Vector w(w_input(0) , w_input(1) , w_input(2) );//陀螺角速度测量值-常漂估计值
		Quaternion A_w(0, w(0), w(1), w(2));
		//状态转移矩阵
		F(0, 1) = w(2), F(0, 2) = -w(1);
		F(1, 0) = -w(2), F(1, 2) = w(0);
		F(2, 0) = w(1), F(2, 1) = -w(0);
		Phi = I + F*dt;
		//姿态估计值
		Quaternion q_;
		//q_ = z + A_w*z*dt*0.5;  //bo1
		Quaternion delta_q_;
		delta_q_ = (A_w*z*dt*0.5);
		delta_q_.w = 1;
		q_ = z*delta_q_;
		normalize(q_);
		
		ekf.TimeUpdate(dt, tr_x0, Phi, Q, w);
		Vector x_est = ekf.State();
		
		//低频噪声
		/*
		Matrix delta_k(3, 6);
		delta_k.M[0][0] = x_est(6), delta_k.M[0][1] = x_est(7),
			delta_k.M[1][2] = x_est(8), delta_k.M[1][3] = x_est(9),
			delta_k.M[2][4] = x_est(10), delta_k.M[2][5] = x_est(11);
		Vector lambda(6);
		lambda(0) = cos(W*dt), lambda(1) = sin(W*dt), lambda(2) = cos(W*dt), lambda(3) = sin(W*dt), lambda(4) = cos(W*dt), lambda(5) = sin(W*dt);
		Vector delta_kAndlambda(3);
		delta_kAndlambda = delta_k*lambda;
		Quaternion delta_kAndlambda2q(1, delta_kAndlambda(0), delta_kAndlambda(1), delta_kAndlambda(2));
		*/
		//q_ = q_*v2q*delta_kAndlambda2q; //bo2
		//q_ = q_*v2q;
		//normalize(q_);
		//量测矩阵
		Matrix H(3, n);
		H(0, 0) = 1, H(1, 1) = 1, H(2, 2) = 1;
		//cout << H << endl;
		/*
		int k = (n - 6) / 6; //傅里叶级数阶数
		for (int i = 1; i <= k; i++)
		{
			H(0, i * 6) = cos(i*W*dt), H(0, i * 6 + 1) = sin(i*W*dt),
				H(1, i * 6 + 2) = cos(i*W*dt), H(1, i * 6 + 3) = sin(i*W*dt),
				H(2, i * 6 + 4) = cos(i*W*dt), H(2, i * 6 + 5) = sin(i*W*dt);
		}
		*/
		//Quaternion g_(x_est[0], x_est[1], x_est[2], x_est[3]); //姿态误差估计值
		g_ = H*x_est;   //Z_k
		ekf.MeasUpdate(g, q_, R, H);
		Vector x_result = ekf.State();
		Matrix P_result = ekf.Cov();
		Quaternion q_result, q_est(1, x_result[0], x_result[1], x_result[2]); //量测更新后的姿态误差
		//normalize(q_est);
		//Quaternion q_0(g_(0), g_(1), g_(2), g_(3));  //估计姿态误差四元数（原来是向量形式）
		//低频噪声
		/*
		Matrix delta_k2(3, 6);
		delta_k2.M[0][0] = x_result(6), delta_k2.M[0][1] = x_result(7),
			delta_k2.M[1][2] = x_result(8), delta_k2.M[1][3] = x_result(9),
			delta_k2.M[2][4] = x_result(10), delta_k2.M[2][5] = x_result(11);
		Vector delta_kAndlambda2(3);
		delta_kAndlambda = delta_k2*lambda;
		Quaternion delta_kAndlambda22q(1, delta_kAndlambda2(0), delta_kAndlambda2(1), delta_kAndlambda2(2));

		q_result = q_*q_est*delta_kAndlambda22q;
		*/
		q_result = q_*q_est;
		normalize(q_result);
		//陀螺常漂校正
		//Vector b_result;
		//Vector delta_b(x_result(3), x_result(4), x_result(5));
		//b_result = b + delta_b;
		e_after = quaternionToEuler(q_result);
		e_before = quaternionToEuler(g);
		//cout << "b估计值x：" << x_result[3] << ",b估计值y：" << x_result[4] << ",b估计值z：" << x_result[5]  << endl;
		//保存估计值姿态
		//e_q_ = quaternionToEuler(q_);
		prev_tm = curr_tm;   //更新上一次的时间
		prev_ms = curr_ms;    //更新上一次的毫秒
		z = q_result; //初始姿态更新
		P0 = P_result; //协方差矩阵更新
		//b = b_result;//陀螺常漂更新
		//b(0) = 0.0005, b(1) = -0.061, b(0) = 0.0002 ;//陀螺常漂更新
		//x_result[0] = q_est.x, x_result[1] = q_est.y, x_result[2] = q_est.z;
		tr_x0 = x_result; //状态量更新

	}
	std::cout << "完成写入！！" << endl;
	std::system("pause");
}
#endif
#endif
#if 0
	//反向滤波初值设定
	Matrix P_b = covarianceMatrixHistory_f.back();
	Vector X_b = stateVectorHistory_f.back();
	EKF ekf_b(n);
	for (int i = data.size() - 1; i >= 0; i--)
	{
		const SensorData& entry = data[i];
		//cout << entry << endl;
		double dt = 0;
		//
		
		//
		dt = timeIntervals[i];
		//cout << dt << endl;
		cout << "----------------------" << endl;
	}
	//cout << P_result << endl;
	


}
#endif