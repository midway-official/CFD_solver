#ifndef FVM_H
#define FVM_H
#include "CFDmath.h"

using namespace std;
//边界条件结构体
struct BoundaryCondition {
    int zoneid;
    int bctype;
    std::string regionName;
};
// 将十六进制的bctype转换为边界类型的描述
std::string getBoundaryType(int bctype);

// 解析网格文件中的边界条件并返回BoundaryCondition的向量
std::vector<BoundaryCondition> parseBoundaryConditions(const std::string& filename);

// 函数声明 
//梯度计算
Field calculateGradient(Field& P, Mesh& mesh);
//散度计算
Field calculateDivergence(Field& U, Mesh& mesh);
#endif // 
