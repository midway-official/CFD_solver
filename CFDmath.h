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
class Point {
public:
    Point(const std::vector<double>& coordinates);

    // 运算符重载
    Point operator+(const Point& other) const;
    Point operator-(const Point& other) const;
    Point operator*(double scalar) const;
    Point operator/(double scalar) const;
    Point operator^(const Point& other) const;
    Point cross(const Point& other) const;
   
    // 计算模长和单位化
    double magnitude() const;
    Point normalize() const;

    // 打印点坐标
    void print() const;

private:
    std::vector<double> coordinates;
};

class Field {
public:
    // 构造函数
    Field(const std::vector<Point>& points);
    Field(const std::vector<Scalar>& scalars);

    // 运算符重载
    Field operator+(const Field& other) const;
    Field operator-(const Field& other) const;
    Field operator*(double scalar) const;
    Field operator/(double scalar) const;

    // 点乘和叉乘
    Scalar dot(const Field& other) const;
    Field cross(const Field& other) const;

    // 标准化和模长
    Field normalize() const;
    Field magnitudes() const; // 返回标量场

    // 打印字段内容
    void print() const;

private:
    std::vector<Point> points;
    std::vector<Scalar> scalars;
    bool isPointField;

    // 通用操作方法
    template<typename Func, typename ResultType = Field>
    ResultType operate(const Field& other, Func func, bool returnScalar = false) const;
};
class Faces {
public:
    // 构造函数
    Faces(const std::vector<std::vector<double>>& nodeCoordinates,
          const std::vector<std::vector<int>>& faceInfo);

    // 打印面的信息
    void printFaces() const;
    // 计算所有面的中心点
    std::vector<Point> calculateAllCenters() const;


private:
    std::vector<std::vector<double>> nodes; // 储存每个节点的三个坐标
    std::vector<std::vector<int>> faces;    // 储存面的信息
};

#endif // CFDMATH_H
