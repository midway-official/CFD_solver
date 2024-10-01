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

    Field dot(const Field& other) const; // 矢量场之间的点乘
    Field cross(const Field& other) const; // 矢量场之间的叉乘

    Field magnitude() const;
    Field normalize() const;

    void print() const;

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

    // 打印面的信息
    void print() const;
    // 计算所有面的中心点
    Field calculateAllCenters() const;
    // 计算所有面法向量 n0 n1 n2右手螺旋确定 由c1指向c0
    Field calculateAllNormals() const;
    //计算所有面面积 返回标量场
    Field calculateAllAreas() const;

private:
    std::vector<std::vector<double>> nodes; // 储存每个节点的三个坐标
    std::vector<std::vector<int>> faces;    // 储存面的信息
};


#endif // CFDMATH_H
