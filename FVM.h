#ifndef FVM_H
#define FVM_H
#include "CFDmath.h"

using namespace std;

struct ZoneInfo {
    int bctype;  // 边界条件类型
    string bctypeText;  // 边界条件类型文本
    bool allPointsMaxX = false;
    bool allPointsMinX = false;
    bool allPointsMaxY = false;
    bool allPointsMinY = false;
    bool allPointsMaxZ = false;
    bool allPointsMinZ = false;
};



class MeshAnalyzer {
public:
    map<int, ZoneInfo> analyzeMesh(const Mesh& mesh);
    void printZoneInfo(const map<int, ZoneInfo>& zoneInfoMap);
private:
    void calculateGlobalExtrema(const Mesh& mesh, double& globalMaxX, double& globalMinX, double& globalMaxY, double& globalMinY, double& globalMaxZ, double& globalMinZ);
};

// 高斯-赛德尔求解器
class GaussSeidel {
public:
    // 构造函数
    GaussSeidel(double tolerance = 0.01, int maxIterations = 1000);

    // 求解函数
    std::vector<double> solve(const std::vector<std::vector<double>>& A, const std::vector<double>& b);

private:
    double tolerance; // 收敛容忍度
    int maxIterations; // 最大迭代次数
};


class GaussSeidel3D {
public:
    // 构造函数
    GaussSeidel3D(double tolerance=10e-8, int maxIterations =1000, double relaxationFactor = 1);

    // 求解函数，A 是系数矩阵，b 是右边的向量矩阵，返回求解结果
    std::vector<std::vector<double>> solve(const std::vector<std::vector<double>>& A,
                                           const std::vector<std::vector<double>>& b);

private:
    double tolerance;         // 容许误差
    int maxIterations;        // 最大迭代次数
    double relaxationFactor;  // 松弛因子
};
class Jacobi3D {
public:
    Jacobi3D(double tolerance = 1e-10, int maxIterations = 1000);
    std::vector<std::vector<double>> solve(const std::vector<std::vector<double>>& A,
                                           const std::vector<double>& b0X,
                                           const std::vector<double>& b0Y,
                                           const std::vector<double>& b0Z);

private:
    double tolerance;
    int maxIterations;
};
// 函数声明 
//组装速度场
std::vector<std::vector<Scalar>> combineVectors(const std::vector<double>& vec1, 
                                                const std::vector<double>& vec2, 
                                                const std::vector<double>& vec3);

//梯度计算(zoneid, bctype, n0, n1, n2, n3, c0, c1)
Field calculateGradient(Field& P, Mesh& mesh,Scalar P_farfield);
//散度计算
Field calculateDivergence(Field& U, Mesh& mesh,Point U_wall);
//离散拉普拉斯项
void separateLaplace( Mesh& mesh, Scalar gramma,vector<vector<Scalar>>& A,vector<Scalar>& b);
void separateLaplaceU(Mesh& mesh, Scalar gramma, vector<vector<Scalar>>& A, vector<vector<Scalar>>& b,Point U_WALL) ;
  
//离散对流项
void separateConvective(Mesh& mesh, Scalar rho, vector<vector<Scalar>>& A, vector<vector<Scalar>>& b,Field Uf,Point Vw );

//离散压力泊松方程
void separateLaplaceP( Mesh& mesh, Scalar gramma,vector<vector<Scalar>>& A,vector<Scalar>& b,vector<vector<Scalar>>& AP);

// 声明计算速度修正的函数
Field computeVelocityCorrection(const Field& p, const std::vector<std::vector<double>>& ap, const Field& vol);

#endif // 