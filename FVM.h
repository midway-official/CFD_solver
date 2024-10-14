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


// 函数声明 
//梯度计算(zoneid, bctype, n0, n1, n2, n3, c0, c1)
Field calculateGradient(Field& P, Mesh& mesh,Scalar P_farfield);
//散度计算
Field calculateDivergence(Field& U, Mesh& mesh,Point U_wall);
#endif // 