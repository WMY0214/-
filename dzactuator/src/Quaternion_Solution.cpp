
#include "Quaternion_Solution.h"
#define SAMPLING_FREQ 20.0f // 采样频率
/**************************************
Date: May 31, 2020
Function: 平方根倒数 求四元数用到
***************************************/
float InvSqrt(float number)
{
  volatile long i;
    volatile float x, y;
    volatile const float f = 1.5F;
    x = number * 0.5F;
    y = number;
    i = * (( long * ) &y);
    i = 0x5f375a86 - ( i >> 1 );
    y = * (( float * ) &i);
    y = y * ( f - ( x * y * y ) );

  return y;
}
/**************************************
Date: May 31, 2020
Function: 四元数解算
***************************************/
volatile float twoKp = 1.0f;     // 2 * proportional gain (Kp)
volatile float twoKi = 0.0f;     // 2 * integral gain (Ki)
volatile float q0 = 1.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;          // quaternion of sensor frame relative to auxiliary frame
volatile float integralFBx = 0.0f,  integralFBy = 0.0f, integralFBz = 0.0f; // integral error terms scaled by Ki
// 基于Mahony互补滤波的四元数姿态解算算法（适用于IMU传感器融合）
void Quaternion_Solution(float gx, float gy, float gz, float ax, float ay, float az)
{
  // 变量声明
  float recipNorm;          // 归一化因子
  float halfvx, halfvy, halfvz; // 重力向量估计值（半计算量）
  float halfex, halfey, halfez; // 重力向量误差
  float qa, qb, qc;         // 四元数临时变量

  // 加速度计数据有效性检查（排除零值干扰）
  if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {
    // 加速度计归一化处理（单位重力向量）
    recipNorm = InvSqrt(ax * ax + ay * ay + az * az);
    ax *= recipNorm;
    ay *= recipNorm;
    az *= recipNorm;      
    
    // 计算当前四元数对应的重力方向（方向余弦矩阵第三行）
    halfvx = q1 * q3 - q0 * q2;  // 估计重力X分量
    halfvy = q0 * q1 + q2 * q3;  // 估计重力Y分量
    halfvz = q0 * q0 - 0.5f + q3 * q3; // 估计重力Z分量
    
    // 计算重力向量误差（叉乘反映方向偏差）
    halfex = (ay * halfvz - az * halfvy); // X轴误差
    halfey = (az * halfvx - ax * halfvz); // Y轴误差
    halfez = (ax * halfvy - ay * halfvx); // Z轴误差

    // 积分反馈项（当前配置twoKi=0未启用）
    if(twoKi > 0.0f) {
      // 积分误差累积（用于消除陀螺仪零偏）
      integralFBx += twoKi * halfex * (1.0f / SAMPLING_FREQ);
      integralFBy += twoKi * halfey * (1.0f / SAMPLING_FREQ);
      integralFBz += twoKi * halfez * (1.0f / SAMPLING_FREQ);
      
      // 应用积分补偿
      gx += integralFBx;
      gy += integralFBy;
      gz += integralFBz;
    }
    else {
      // 重置积分项防止饱和（当前配置）
      integralFBx = integralFBy = integralFBz = 0.0f;
    }

    // 比例反馈补偿（当前twoKp=1.0）
    gx += twoKp * halfex;  // X轴角速度补偿
    gy += twoKp * halfey;  // Y轴角速度补偿
    gz += twoKp * halfez;  // Z轴角速度补偿
  }

  // 四元数微分方程积分（20Hz采样周期）
  gx *= (0.5f * (1.0f / SAMPLING_FREQ)); // 转换为角度增量
  gy *= (0.5f * (1.0f / SAMPLING_FREQ));
  gz *= (0.5f * (1.0f / SAMPLING_FREQ));
  
  // 四元数更新（一阶龙格库塔法）
  qa = q0; qb = q1; qc = q2;
  q0 += (-qb * gx - qc * gy - q3 * gz);
  q1 += (qa * gx + qc * gz - q3 * gy);
  q2 += (qa * gy - qb * gz + q3 * gx);
  q3 += (qa * gz + qb * gy - qc * gx); 

  // 四元数归一化（防止数值漂移）
  recipNorm = InvSqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
  q0 *= recipNorm;
  q1 *= recipNorm;
  q2 *= recipNorm;
  q3 *= recipNorm;

  // 输出到IMU数据结构（ROS标准格式）
  Mpu6050.orientation.w = q0;
  Mpu6050.orientation.x = q1;
  Mpu6050.orientation.y = q2;
  Mpu6050.orientation.z = q3;
}
