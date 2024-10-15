#include "CFDmath.h"  // 数学库头文件
#include "FVM.h"
#include <map>
#include <fstream>
#include <iostream>
#include <regex>
#include <vector>

using namespace std;
// 定义边界条件类型映射
map<int, string> boundaryTypes = {
    {2, "interior"},
    {3, "wall"},
    {4, "pressure-inlet"},
    {5, "pressure-outlet"},
    {7, "symmetry"},
    {8, "periodic-shadow"},
    {9, "pressure-far-field"},
    {11, "velocity-inlet"},
    {12, "periodic"},
    {14, "fan, porous-jump, radiator"},
    {20, "mass-flow-inlet"},
    {24, "Interface"},
    {31, "parent"},
    {36, "outflow"},
    {37, "axis"}
};
//组装速度场
std::vector<std::vector<Scalar>> combineVectors(const std::vector<double>& vec1, 
                                                const std::vector<double>& vec2, 
                                                const std::vector<double>& vec3) {
    // 检查输入向量的长度是否相同
    if (vec1.size() != vec2.size() || vec1.size() != vec3.size()) {
        throw std::invalid_argument("All vectors must have the same length.");
    }

    // 创建一个三列的vector<vector<Scalar>>
    std::vector<std::vector<Scalar>> combined(vec1.size(), std::vector<Scalar>(3));

    for (size_t i = 0; i < vec1.size(); ++i) {
        combined[i][0] = vec1[i]; // 第一列
        combined[i][1] = vec2[i]; // 第二列
        combined[i][2] = vec3[i]; // 第三列
    }

    return combined;
}
// 计算整个网格的 x, y, z 最大最小值
void MeshAnalyzer::calculateGlobalExtrema(const Mesh& mesh, double& globalMaxX, double& globalMinX, double& globalMaxY, double& globalMinY, double& globalMaxZ, double& globalMinZ) {
    globalMaxX = globalMaxY = globalMaxZ = -numeric_limits<double>::infinity();
    globalMinX = globalMinY = globalMinZ = numeric_limits<double>::infinity();

    size_t numPoints = mesh.numberOfPoints();
    for (size_t i = 0; i < numPoints; ++i) {
        const vector<double>& point = mesh.getPoint(i);

        globalMaxX = max(globalMaxX, point[0]);
        globalMinX = min(globalMinX, point[0]);
        globalMaxY = max(globalMaxY, point[1]);
        globalMinY = min(globalMinY, point[1]);
        globalMaxZ = max(globalMaxZ, point[2]);
        globalMinZ = min(globalMinZ, point[2]);
    }
}

map<int, ZoneInfo> MeshAnalyzer::analyzeMesh(const Mesh& mesh) {
    map<int, ZoneInfo> zoneInfoMap;

    double globalMaxX, globalMinX, globalMaxY, globalMinY, globalMaxZ, globalMinZ;
    calculateGlobalExtrema(mesh, globalMaxX, globalMinX, globalMaxY, globalMinY, globalMaxZ, globalMinZ);

    size_t numFaces = mesh.numberOfFaces();
    for (size_t i = 0; i < numFaces; ++i) {
        const vector<int>& face = mesh.getFace(i);
        int zoneid = face[0];
        int bctype = face[1];
        int n0 = face[2], n1 = face[3], n2 = face[4], n3 = face[5];

        // 提取该 zone 的 4 个点的坐标
        const vector<double>& p0 = mesh.getPoint(n0 - 1);
        const vector<double>& p1 = mesh.getPoint(n1 - 1);
        const vector<double>& p2 = mesh.getPoint(n2 - 1);
        const vector<double>& p3 = mesh.getPoint(n3 - 1);

        // 初始化或更新该 zone 的信息
        if (zoneInfoMap.find(zoneid) == zoneInfoMap.end()) {
            zoneInfoMap[zoneid] = ZoneInfo{bctype, boundaryTypes[bctype]};
        }
        ZoneInfo& zone = zoneInfoMap[zoneid];

        // 检查是否所有点的 x 坐标相同，且等于全局最大值或最小值
        bool sameX = (p0[0] == p1[0] && p0[0] == p2[0] && p0[0] == p3[0]);
        bool sameY = (p0[1] == p1[1] && p0[1] == p2[1] && p0[1] == p3[1]);
        bool sameZ = (p0[2] == p1[2] && p0[2] == p2[2] && p0[2] == p3[2]);

        if (sameX && p0[0] == globalMaxX) {
            zone.allPointsMaxX = true;
        }
        if (sameX && p0[0] == globalMinX) {
            zone.allPointsMinX = true;
        }
        if (sameY && p0[1] == globalMaxY) {
            zone.allPointsMaxY = true;
        }
        if (sameY && p0[1] == globalMinY) {
            zone.allPointsMinY = true;
        }
        if (sameZ && p0[2] == globalMaxZ) {
            zone.allPointsMaxZ = true;
        }
        if (sameZ && p0[2] == globalMinZ) {
            zone.allPointsMinZ = true;
        }
    }

    return zoneInfoMap;
}

void MeshAnalyzer::printZoneInfo(const map<int, ZoneInfo>& zoneInfoMap) {
    for (const auto& [zoneid, info] : zoneInfoMap) {
        cout << "Zone ID: " << zoneid << ", BCType: " << info.bctype << " (" << info.bctypeText << ") - ";

        if (info.allPointsMaxX) {
            cout << "Max X ";
        } else if (info.allPointsMinX) {
            cout << "Min X ";
        }
        if (info.allPointsMaxY) {
            cout << "Max Y ";
        } else if (info.allPointsMinY) {
            cout << "Min Y ";
        }
        if (info.allPointsMaxZ) {
            cout << "Max Z ";
        } else if (info.allPointsMinZ) {
            cout << "Min Z ";
        }

        cout << endl;
    }
}

// 构造函数
GaussSeidel::GaussSeidel(double tolerance, int maxIterations)
    : tolerance(tolerance), maxIterations(maxIterations) {}

std::vector<double> GaussSeidel::solve(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
    int n = A.size(); // 系数矩阵的行数
    std::vector<double> x(n, 0.0); // 初始化解向量

    for (int k = 0; k < maxIterations; ++k) {
        std::vector<double> x_old = x; // 记录上一次的解

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;

            // 计算 A[i][j] * x[j] 的和
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum += A[i][j] * x[j];
                }
            }

            // 更新当前变量
            x[i] = (b[i] - sum) / A[i][i];
        }

        // 计算残差 r = b - A * x
        double residual_sum = 0.0;
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                sum += A[i][j] * x[j];
            }
            residual_sum += fabs(b[i] - sum); // 计算残差的绝对值，并累加
        }

        // 输出当前迭代的残差总和
        //std::cout << "Iteration " << k + 1 << ", Residual sum: " << residual_sum << std::endl;

        // 检查收敛性：比较新旧解之间的变化量
        double norm = 0.0;
        for (int i = 0; i < n; ++i) {
            norm += pow(x[i] - x_old[i], 2);
        }
        norm = sqrt(norm);

        if (norm < tolerance) {
            //std::cout << "Converged after " << k + 1 << " iterations." << std::endl;
            return x;
        }
    }

    std::cout << "Max iterations reached." << std::endl;
    return x; // 如果没有收敛，返回最后的解
}


// 构造函数
GaussSeidel3D::GaussSeidel3D(double tolerance, int maxIterations, double relaxationFactor)
    : tolerance(tolerance), maxIterations(maxIterations), relaxationFactor(relaxationFactor) {}

// 求解函数
std::vector<std::vector<double>> GaussSeidel3D::solve(const std::vector<std::vector<double>>& A,
                                                      const std::vector<std::vector<double>>& b) {
    int n = A.size(); // 系数矩阵的行数
    std::vector<std::vector<double>> x(n, std::vector<double>(3, 0.0)); // 初始化解向量为0

    for (int k = 0; k < maxIterations; ++k) {
        std::vector<std::vector<double>> x_old = x; // 记录上一次的解

        for (int i = 0; i < n; ++i) {
            // 计算 A[i][j] * x[j] 的和
            double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sumX += A[i][j] * x[j][0];
                    sumY += A[i][j] * x[j][1];
                    sumZ += A[i][j] * x[j][2];
                }
            }

            // 更新当前变量，加入松弛因子
            x[i][0] = (1 - relaxationFactor) * x_old[i][0] + 
                       (relaxationFactor / A[i][i]) * (b[i][0] - sumX); // 更新 x 分量
            x[i][1] = (1 - relaxationFactor) * x_old[i][1] + 
                       (relaxationFactor / A[i][i]) * (b[i][1] - sumY); // 更新 y 分量
            x[i][2] = (1 - relaxationFactor) * x_old[i][2] + 
                       (relaxationFactor / A[i][i]) * (b[i][2] - sumZ); // 更新 z 分量
        }

        // 计算残差 r = b - A * x
        double residual_sum = 0.0;
        for (int i = 0; i < n; ++i) {
            double AxX = 0.0, AxY = 0.0, AxZ = 0.0;
            for (int j = 0; j < n; ++j) {
                AxX += A[i][j] * x[j][0];
                AxY += A[i][j] * x[j][1];
                AxZ += A[i][j] * x[j][2];
            }
            residual_sum += fabs(b[i][0] - AxX) + fabs(b[i][1] - AxY) + fabs(b[i][2] - AxZ);
        }

        // 输出当前迭代的残差总和
        //std::cout << "Iteration " << k + 1 << ", Residual sum: " << residual_sum << std::endl;

        // 检查收敛性：比较新旧解之间的变化量
        double norm = 0.0;
        for (int i = 0; i < n; ++i) {
            norm += pow(x[i][0] - x_old[i][0], 2) + pow(x[i][1] - x_old[i][1], 2) + pow(x[i][2] - x_old[i][2], 2);
        }
        norm = sqrt(norm);

        if (norm < tolerance) {
            //std::cout << "Converged after " << k + 1 << " iterations." << std::endl;
            return x; // 返回解向量
        }
    }

    std::cout << "Max iterations reached." << std::endl;
    return x; // 如果没有收敛，返回最后的解
}

Jacobi3D::Jacobi3D(double tolerance, int maxIterations)
    : tolerance(tolerance), maxIterations(maxIterations) {}

std::vector<std::vector<double>> Jacobi3D::solve(const std::vector<std::vector<double>>& A,
                                                  const std::vector<double>& b0X,
                                                  const std::vector<double>& b0Y,
                                                  const std::vector<double>& b0Z) {
    int n = A.size(); // 系数矩阵的行数
    std::vector<std::vector<double>> x(n, std::vector<double>(3, 0.0)); // 初始化解向量为0

    for (int k = 0; k < maxIterations; ++k) {
        std::vector<std::vector<double>> x_new = x; // 新的解向量

        for (int i = 0; i < n; ++i) {
            double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sumX += A[i][j] * x[j][0]; // 使用前一次的解
                    sumY += A[i][j] * x[j][1];
                    sumZ += A[i][j] * x[j][2];
                }
            }

            // 更新当前变量
            x_new[i][0] = (b0X[i] - sumX) / A[i][i]; // 更新 x 分量
            x_new[i][1] = (b0Y[i] - sumY) / A[i][i]; // 更新 y 分量
            x_new[i][2] = (b0Z[i] - sumZ) / A[i][i]; // 更新 z 分量
        }

        // 计算残差 r = b - A * x
        double residual_sum = 0.0;
        for (int i = 0; i < n; ++i) {
            double AxX = 0.0, AxY = 0.0, AxZ = 0.0;
            for (int j = 0; j < n; ++j) {
                AxX += A[i][j] * x_new[j][0];
                AxY += A[i][j] * x_new[j][1];
                AxZ += A[i][j] * x_new[j][2];
            }
            residual_sum += fabs(b0X[i] - AxX) + fabs(b0Y[i] - AxY) + fabs(b0Z[i] - AxZ);
        }

        // 输出当前迭代的残差总和
        std::cout << "Iteration " << k + 1 << ", Residual sum: " << residual_sum << std::endl;

        // 检查收敛性：比较新旧解之间的变化量
        double norm = 0.0;
        for (int i = 0; i < n; ++i) {
            norm += pow(x_new[i][0] - x[i][0], 2) + pow(x_new[i][1] - x[i][1], 2) + pow(x_new[i][2] - x[i][2], 2);
        }
        norm = sqrt(norm);

        if (norm < tolerance) {
            std::cout << "Converged after " << k + 1 << " iterations." << std::endl;
            return x_new; // 返回新的解向量
        }

        x = x_new; // 更新解向量
    }

    std::cout << "Max iterations reached." << std::endl;
    return x; // 如果没有收敛，返回最后的解
}
//函数区域
// 函数：计算压力梯度
Field calculateGradient(Field& P, Mesh& mesh,Scalar P_farfield) {
    Field gradientP(mesh.numberOfCells(), vector<double>{0, 0, 0});
    Field facearea = mesh.calculateAllfaceAreas();
    Field facenorm = mesh.calculateAllfaceNormals();
    Field cellVol = mesh.calculateAllcellVolumes();

    for (size_t i = 0; i < mesh.numberOfFaces(); ++i) {
        int bctype = mesh.getFace(i)[1]; // 获取bctype

        switch(bctype) {
            case 2: { // 内部面
                size_t c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                size_t c1 = mesh.getFace(i)[7] - 1;//(zoneid, bctype, n0, n1, n2, n3, c0, c1)

                // 确保 c0 和 c1 在有效范围内
                if (c0 >= P.size() || c1 >= P.size()) {
                    cerr << "Error: Cell index out of range. c0: " << c0 << ", c1: " << c1 << ", P size: " << P.size() << endl;
                    throw out_of_range("Cell超出索引.");
                }

                Scalar Pf = -(P.scalarAt(c0) + P.scalarAt(c1)) / 2.0; // 法向量指向c0，故加负号表示朝外
                double A = facearea.scalarAt(i); // 当前面面积

                Point gradientPc0 = gradientP.pointAt(c0);
                Point gradientPc1 = gradientP.pointAt(c1);

                gradientP.setPoint(c0, gradientPc0 + facenorm.pointAt(i) * Pf * A);
                gradientP.setPoint(c1, gradientPc1 + facenorm.pointAt(i) * Pf * A * (-1.0)); // c0 c1外向面法向量相反
                break;
            }
            case 9: { // 压力远场或描述压强
                size_t c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
               //c1为0(zoneid, bctype, n0, n1, n2, n3, c0, c1)

                // 确保 c0 和 c1 在有效范围内
                if (c0 >= P.size() ) {
                    cerr << "Error: Cell index out of range. c0: " << c0 << ", c1: " << "压力远场" << ", P size: " << P.size() << endl;
                    throw out_of_range("Cell超出索引.");
                }

                Scalar Pf = -(P.scalarAt(c0) + P_farfield) / 2.0; // 法向量指向c0，故加负号表示朝外
                double A = facearea.scalarAt(i); // 当前面面积

                Point gradientPc0 = gradientP.pointAt(c0);
                

                gradientP.setPoint(c0, gradientPc0 + facenorm.pointAt(i) * Pf * A);
                
                break;
            }
            case 36: { // 零压梯度 出流条件
                size_t c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
               //c1为0(zoneid, bctype, n0, n1, n2, n3, c0, c1)

                // 确保 c0 和 c1 在有效范围内
                if (c0 >= P.size() ) {
                    cerr << "Error: Cell index out of range. c0: " << c0 << ", c1: " << "零压梯度" << ", P size: " << P.size() << endl;
                    throw out_of_range("Cell超出索引.");
                }

                Scalar Pf = -P.scalarAt(c0); // 法向量指向c0，故加负号表示朝外
                double A = facearea.scalarAt(i); // 当前面面积

                Point gradientPc0 = gradientP.pointAt(c0);
                

                gradientP.setPoint(c0, gradientPc0 + facenorm.pointAt(i) * Pf * A);
                
                break;
            }
            // 其他面类型可以在此添加不同的处理逻辑
            default:
                // 对于其他类型的面，跳过
                //wall不用处理
                continue;
        }
    }


    return gradientP / cellVol;
}

// 函数：计算散度
Field calculateDivergence(Field& U, Mesh& mesh, Point U_wall) {
    Field divergenceU(mesh.numberOfCells(), 0.0);
    Field facearea = mesh.calculateAllfaceAreas();
    Field facenorm = mesh.calculateAllfaceNormals();
    Field cellVol = mesh.calculateAllcellVolumes();

    for (size_t i = 0; i < mesh.numberOfFaces(); ++i) {
        int bctype = mesh.getFace(i)[1]; // 获取bctype
        switch(bctype) {
            case 2:{ // 内部面
                size_t c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                size_t c1 = mesh.getFace(i)[7] - 1;

                Point Uf = (U.pointAt(c0) + U.pointAt(c1)) * (-0.5); // 计算面上的速度，加面法向量指入c0负号

                Scalar A = facearea.scalarAt(i); // 当前面面积

                // 散度更新
                Scalar divergenceContributionC0 = A * (Uf.dot(facenorm.pointAt(i)));
                Scalar divergenceContributionC1 = A * (Uf.dot(facenorm.pointAt(i)));

                divergenceU.setScalar(c0, divergenceU.scalarAt(c0) + divergenceContributionC0);
                divergenceU.setScalar(c1, divergenceU.scalarAt(c1) - divergenceContributionC1); // 法向量相反
                break;
            }
            case 9:{ // 压力远场，描述压强
                size_t c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                
                Point p0(vector<Scalar>{0.0, 0.0, 0.0});
                Point Uf = (U.pointAt(c0) + p0) * (-0.5); // 计算面上的速度，加面法向量指入c0负号

                Scalar A = facearea.scalarAt(i); // 当前面面积

                // 散度更新
                Scalar divergenceContributionC0 = A * (Uf.dot(facenorm.pointAt(i)));
                divergenceU.setScalar(c0, divergenceU.scalarAt(c0) + divergenceContributionC0);
                break;
            }
            case 36:{ // 零压梯度 出流条件
                size_t c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                
                Point p0(vector<Scalar>{0.0, 0.0, 0.0});
                Point Uf = (U.pointAt(c0) + p0) * (-0.5); // 计算面上的速度，加面法向量指入c0负号

                Scalar A = facearea.scalarAt(i); // 当前面面积

                // 散度更新
                Scalar divergenceContributionC0 = A * (Uf.dot(facenorm.pointAt(i)));
                divergenceU.setScalar(c0, divergenceU.scalarAt(c0) + divergenceContributionC0);
                break;
            }
            case 3:{ // wall壁面
                size_t c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                
                Point Uf = (U.pointAt(c0) + U_wall) * (-0.5); // 计算面上的速度，加面法向量指入c0负号

                Scalar A = facearea.scalarAt(i); // 当前面面积

                // 散度更新
                Scalar divergenceContributionC0 = A * (Uf.dot(facenorm.pointAt(i)));
                divergenceU.setScalar(c0, divergenceU.scalarAt(c0) + divergenceContributionC0);
                break;
            }
            default:
                // 对于其他类型的面，跳过
                continue;
        }
    }



    // 确保在for循环结束后返回结果
    return divergenceU / cellVol;
}
//使用时右侧乘v
void separateLaplace(Mesh& mesh, Scalar gramma, vector<vector<Scalar>>& A, vector<Scalar>& b) {
    int n = mesh.numberOfCells();
    
    Field facearea = mesh.calculateAllfaceAreas();
    Field facenorm = mesh.calculateAllfaceNormals();
    Field cellVol = mesh.calculateAllcellVolumes();
    Field cellCenter = mesh.calculateAllcellCenters();
    Field faceCenter = mesh.calculateAllfaceCenters();

    
    for (size_t i = 0; i < mesh.numberOfFaces(); i++) {
        int bctype = mesh.getFace(i)[1]; // 获取bctype
   
        
        // 在 switch 语句外部声明变量
        size_t c0, c1; // 定义 c0 和 c1
        Scalar distance, area, Vol, dS;

        switch (bctype) {
            case 2: // 内部面
                c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                c1 = mesh.getFace(i)[7] - 1;

                distance = (cellCenter.pointAt(c0) - cellCenter.pointAt(c1)).magnitude();
                area = facearea.scalarAt(i);
                
                dS = -gramma * area / (distance );

                A[c0][c0] += dS;
                A[c1][c0] -= dS;
                A[c0][c1] -= dS;
                A[c1][c1] += dS;
                break;

            case 9: // 压力远场
                c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                
                distance = (cellCenter.pointAt(c0) - faceCenter.pointAt(i)).magnitude();
                area = facearea.scalarAt(i);
              
                dS = -gramma * area / (distance );

                A[c0][c0] += dS;
                b[c0] += 0 * dS; // 边界上压强为0
                break;

            case 3: // 壁面
            case 36: // 零压梯度
                c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
               
                A[c0][c0] += 0 * dS; // 不影响
                b[c0] += 0 * dS;
                break;

            default:
                break;
        }
    }
}


void separateLaplaceU(Mesh& mesh, Scalar gramma, vector<vector<Scalar>>& A, vector<vector<Scalar>>& b,Point U_WALL) {
    int n = mesh.numberOfCells();
    
    Field facearea = mesh.calculateAllfaceAreas();
    Field facenorm = mesh.calculateAllfaceNormals();
    Field cellVol = mesh.calculateAllcellVolumes();
    Field cellCenter = mesh.calculateAllcellCenters();
    Field faceCenter = mesh.calculateAllfaceCenters();

    
    for (size_t i = 0; i < mesh.numberOfFaces(); i++) {
        int bctype = mesh.getFace(i)[1]; // 获取bctype
   
        
        // 在 switch 语句外部声明变量
        size_t c0, c1; // 定义 c0 和 c1
        Scalar distance, area, Vol, dS;

        switch (bctype) {
            case 2: // 内部面
                c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                c1 = mesh.getFace(i)[7] - 1;

                distance = (cellCenter.pointAt(c0) - cellCenter.pointAt(c1)).magnitude();
                area = facearea.scalarAt(i);
                
                dS = -gramma * area / (distance );

                A[c0][c0] += dS;
                A[c1][c0] -= dS;
                A[c0][c1] -= dS;
                A[c1][c1] += dS;
                break;

            case 9: // 压力远场
            case 36: // 零压梯度
                c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                
                distance = (cellCenter.pointAt(c0) - faceCenter.pointAt(i)).magnitude();
                area = facearea.scalarAt(i);
              
                dS = -gramma * area / (distance );

                A[c0][c0] += dS;
                b[c0][0] += 0 * dS; // 边界上压强为0
                b[c0][1] += 0 * dS; // 边界上压强为0
                b[c0][2] += 0 * dS; // 边界上压强为0
                break;

            case 3: // 壁面
            
                c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
               
                A[c0][c0] += dS;
                b[c0][0] += U_WALL.getCoordinate(0) * dS; // 边界上压强为0
                b[c0][1] +=  U_WALL.getCoordinate(1) * dS;
                b[c0][2] +=  U_WALL.getCoordinate(2) * dS;
                break;

            default:
                break;
        }
    }
}

//使用时右侧乘v
void separateConvective(Mesh& mesh, Scalar rho, vector<vector<Scalar>>& A, vector<vector<Scalar>>& b,Field Uf,Point Vw ){
     int n = mesh.numberOfCells();
    
    Field facearea = mesh.calculateAllfaceAreas();
    Field facenorm = mesh.calculateAllfaceNormals();
    Field cellVol = mesh.calculateAllcellVolumes();
    Field cellCenter = mesh.calculateAllcellCenters();
    Field faceCenter = mesh.calculateAllfaceCenters();

   
    for (size_t i = 0; i < mesh.numberOfFaces(); i++) {
        int bctype = mesh.getFace(i)[1]; // 获取bctype
   
        
        // 在 switch 语句外部声明变量
        size_t c0, c1; // 定义 c0 和 c1
        Scalar distance, area, Vol, dS,mf;
        Point u0,u1;
        switch (bctype) {
            case 2: // 内部面
                c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                c1 = mesh.getFace(i)[7] - 1;
                u0=  Uf.pointAt(c0);
                u1=  Uf.pointAt(c1);
                area = facearea.scalarAt(i);
                if ((facenorm.pointAt(i)).dot(u1) >= 0)//速度从c1流向c0
                {
                    mf=-u1.dot(facenorm.pointAt(i))*area*rho;//流入为负
                    }else{
                    mf=-u0.dot(facenorm.pointAt(i))*area*rho;//流出为正
                }
                if (mf >=0)
                {
                    A[c0][c0] +=mf;
                    A[c1][c0] -=mf;
                }
                else
                {
                    A[c0][c1] +=mf;
                    A[c1][c1] -=mf;
                }
                
               
                
                break;
            case 36: //零压梯度
            case 9: // 压力远场
                c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                
                u0=  Uf.pointAt(c0);
                u1=  Point(vector<Scalar>(3,0));
                area = facearea.scalarAt(i);
                if ((facenorm.pointAt(i)).dot(u0) >= 0)//速度从c1流向c0
                {
                    mf=-u1.dot(facenorm.pointAt(i))*area*rho;//流入为负
                    }else{
                    mf=-u0.dot(facenorm.pointAt(i))*area*rho;//流出为正
                }
                if (mf >=0)
                {
                    A[c0][c0] +=mf;
                    
                }
                else
                {
                    
                    b[c0][0]-=mf*u1.getCoordinate(0);
                    b[c0][1]-=mf*u1.getCoordinate(1);
                    b[c0][2]-=mf*u1.getCoordinate(2);
                }
                

            case 3: // 壁面
           
                 c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                
                u0=  Uf.pointAt(c0);
                u1=  Vw;
                area = facearea.scalarAt(i);
                if ((facenorm.pointAt(i)).dot(u0) >= 0)//速度从c1流向c0
                {
                    mf=-u1.dot(facenorm.pointAt(i))*area*rho;//流入为负
                    }else{
                    mf=-u0.dot(facenorm.pointAt(i))*area*rho;//流出为正
                }
                if (mf >=0)
                {
                    A[c0][c0] +=mf;
                    
                }
                else
                {
                    
                    b[c0][0]-=mf*u1.getCoordinate(0);
                    b[c0][1]-=mf*u1.getCoordinate(1);
                    b[c0][2]-=mf*u1.getCoordinate(2);
                }
                break;

            default:
                break;
        }
    }
    

}

void separateLaplaceP(Mesh& mesh, Scalar gramma, vector<vector<Scalar>>& A, vector<Scalar>& b,vector<vector<Scalar>>& AP) {
    int n = mesh.numberOfCells();
    
    Field facearea = mesh.calculateAllfaceAreas();
    Field facenorm = mesh.calculateAllfaceNormals();
    Field cellVol = mesh.calculateAllcellVolumes();
    Field cellCenter = mesh.calculateAllcellCenters();
    Field faceCenter = mesh.calculateAllfaceCenters();

    
    for (size_t i = 0; i < mesh.numberOfFaces(); i++) {
        int bctype = mesh.getFace(i)[1]; // 获取bctype
   
        
        // 在 switch 语句外部声明变量
        size_t c0, c1; // 定义 c0 和 c1
        Scalar distance, area, Vol, dS;

        switch (bctype) {
            case 2: // 内部面
                c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                c1 = mesh.getFace(i)[7] - 1;

                distance = (cellCenter.pointAt(c0) - cellCenter.pointAt(c1)).magnitude();
                area = facearea.scalarAt(i);
                
                dS = -gramma * area / (distance );

                A[c0][c0] +=AP[c0][c0]*dS /cellVol.scalarAt(c0) ;
                A[c1][c0] -= AP[c1][c1]*dS /cellVol.scalarAt(c1);
                A[c0][c1] -= AP[c0][c0]*dS /cellVol.scalarAt(c0);
                A[c1][c1] += AP[c1][c1]*dS /cellVol.scalarAt(c1);
                break;

            /*case 9: // 压力远场
                c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
                
                distance = (cellCenter.pointAt(c0) - faceCenter.pointAt(i)).magnitude();
                area = facearea.scalarAt(i);
              
                dS = -gramma * area / (distance );

                A[c0][c0] += AP[c0][c0]*dS /cellVol.scalarAt(c0);
                b[c0] += 0 *AP[c0][c0]*dS /cellVol.scalarAt(c0); // 边界上压强为0
                break;

            case 3: // 壁面
            case 36: // 零压梯度
                c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
               
                A[c0][c0] += 0 * AP[c0][c0]*dS /cellVol.scalarAt(c0); // 不影响
                b[c0] += 0 * AP[c0][c0]*dS /cellVol.scalarAt(c0);
                break;*/

            default:
                break;
        }
    }
}

// 计算速度修正的实现
Field computeVelocityCorrection(const Field& p, const std::vector<std::vector<double>>& ap, const Field& vol) {
    // 初始化一个和 p_ 大小相同的矢量场 velocityCorr，用于存储计算结果
    Field velocityCorr(p.size(), vector<Scalar>{0.0, 0.0, 0.0});

    // 遍历场的所有点
    for (size_t i = 0; i < p.size(); ++i) {
        // 获取第 i 个点的压力修正矢量 p_i
        Point p_i = p.pointAt(i);

        // 获取第 i 个点的体积和 ap 矩阵的主对角线元素 ap[i][i]
        double vol_i = vol.scalarAt(i);
        double ap_ii = ap[i][i];  // 主对角元

        // 计算速度修正：(-vol_i / ap_ii) * p_i
        double scale = -vol_i / ap_ii;
        Point velocityCorr_i = p_i * scale;

        // 更新 velocityCorr 场中的第 i 个点
        velocityCorr.setPoint(i, velocityCorr_i);
    }

    // 返回计算好的速度修正场
    return velocityCorr;
}