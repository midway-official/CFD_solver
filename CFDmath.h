#ifndef CFDMATH_H
#define CFDMATH_H
#include <vector>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <regex>
#include <cmath>
#include <unordered_set>
using namespace std;
using std::ifstream;   // 简化使用 ifstream
using std::stringstream; // 简化使用 stringstream
using std::string;     // 简化使用 string
using std::vector;     // 简化使用 vector
// 假设 Scalar 是 double 的别名
using Scalar = double;
//定义读取msh文件nodes信息函数
void parseNodes(const std::string &filename, std::vector<std::vector<double>> &coordinates);
//定义读取msh文件faces信息函数
void parseFaces(const std::string &filename, std::vector<std::vector<int>> &faces);
//定义读取msh文件facetree信息函数
void parseFacetree(const std::string &filename, std::vector<std::vector<int>> &facetree);
//定义生成celldata信息
std::vector<std::vector<int>> processCelldata(std::vector<std::vector<int>> faceData, const std::vector<std::vector<int>>& facetree);
// 打印二维整型向量中的最大值
void printMaxValue(const std::vector<std::vector<int>>& vec);

// 打印节点数据
void printNodeData(const std::vector<std::vector<double>>& nodeData);

// 打印面数据
void printData(const std::vector<std::vector<int>>& faceData);

// 打印面集合
void printFaces(const std::vector<std::vector<int>>& faces);

// 打印子向量的最大长度
void printMaxSubvectorLength(const std::vector<std::vector<int>>& vec);

//点类定义
class Point {
public:
    Point(const std::vector<double>& coordinates);

    // 运算符重载
    Point operator+(const Point& other) const;
    Point operator-(const Point& other) const;
    Point operator*(double scalar) const;
    Point operator/(double scalar) const;
  //叉乘^
    Point operator^(const Point& other) const;
//叉乘cross
    Point cross(const Point& other) const;
//点乘
   Scalar dot(const Point& other) const;
    // 计算模长和单位化
    double magnitude() const;
    Point normalize() const;

    // 打印点坐标
    void print() const;

private:
    std::vector<double> coordinates;
};
//场 矢量场或标量场
class Field {
public:
    Field(const std::vector<Point>& points);
    Field(const std::vector<Scalar>& scalars);

    Field operator+(const Field& other) const;
    Field operator-(const Field& other) const;
    Field operator*(double scalar) const; // 矢量场与标量的乘法
    Field operator*(const Field& other) const; // 标量场与矢量场的乘法
    //友元函数，标量乘场
    friend Field operator*(double scalar, const Field& field);

    Field dot(const Field& other) const; // 矢量场之间的点乘
    Field cross(const Field& other) const; // 矢量场之间的叉乘

    Field magnitude() const;//求模
    Field normalize() const;//标准化
    void print() const ;//打印场
   // 针对点的运算符重载 索引从0开始
    Point& pointAt(size_t index);
    const Point& pointAt(size_t index) const;

    // 针对标量的运算符重载 索引
    Scalar& scalarAt(size_t index);
    const Scalar& scalarAt(size_t index) const;

private:
    void checkCompatibility(const Field& other) const;

    bool isPointField;
    std::vector<Point> points;
    std::vector<Scalar> scalars;
};
//面集faces
class Faces {
public:
    // 构造函数
    Faces(const std::vector<std::vector<double>>& nodeCoordinates,
          const std::vector<std::vector<int>>& faceInfo);
     Faces();
    // 打印面的信息
    void print() const;
    // 计算所有面的中心点返回矢量场
    Field calculateAllCenters() const;
    // 计算所有面法向量 n0 n1 n2右手螺旋确定 由c1指向c0返回矢量场
    Field calculateAllNormals() const;
    //计算所有面面积 返回标量场
    Field calculateAllAreas() const;
    // 访问索引为 i 的面
    const std::vector<int>& getFace(size_t i) const;

    // 修改 faces 中第一个元素为 i 的行的第二个元素
    void setBctypeForFace(int zoneIndex, int bctype);

private:
    std::vector<std::vector<double>> nodes; // 储存每个节点的三个坐标
    std::vector<std::vector<int>> faces;    // 储存面的信息
};
//cells单位体集
class Cells {
public:
     Cells();
     Cells(const std::vector<std::vector<double>>& nodeCoordinates,
        const std::vector<std::vector<int>>& cellsdata);
    //计算所有单元cell中心点
     Field calculateAllCenters() const;
    //计算所有单元cell体积
     Field calculateAllVolumes() const;
private:
    std::vector<std::vector<double>> nodes; // 储存每个节点的三个坐标
    std::vector<std::vector<int>> cells;    // 储存cell单元体的信息
};
//网格对象Mesh
class Mesh {
public:
    //构造函数
    Mesh(const std::string &filename);
    //计算所有面中心
    Field calculateAllfaceCenters() const;
    //计算所有面法向量
    Field calculateAllfaceNormals() const;
    //计算所有面面积
    Field calculateAllfaceAreas() const;
    //计算所有体cell体心
    Field calculateAllcellCenters() const;
    //计算所有体体积
    Field calculateAllcellVolumes() const;
    //获得第i个面
    const std::vector<int>& getFace(size_t i) const;
    //设置zoneid的边界条件
    void setBctypeForFace(int zoneIndex, int bctype);

private:
    std::vector<std::vector<double>> nodes;
    Faces faces;
    Cells cells;
};

#endif // CFDMATH_H
