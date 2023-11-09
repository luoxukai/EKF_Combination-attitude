# Kalman-filtering
The EKF is used to combine the observation data of the star sensor and the measurement data of the gyroscope
样例数据为一颗卫星上双星敏所观测的四元数、双陀螺测量角速度、以及陀螺温度

主函数状态量为姿态误差，由公式可推导姿态误差可降阶为三项，取矢量部分。另代码状态量部分可扩展加入陀螺常漂误差、星敏低频噪声系数误差，分别为3维和6n维（n为低频噪声所取傅里叶级数阶数）

另可加入反向滤波部分进行融合滤波
