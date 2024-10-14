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