#ifndef CFDMATH_H  // 防止头文件重复包含
#define CFDMATH_H

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <regex>
#include <cmath>
// 假设 Scalar 是 double 的别名
using Scalar = double;

// 打印二维数组大小
template <typename T>
inline void print2DArraySize(const std::vector<std::vector<T>>& arr) {
    std::cout << "Size of array: " << arr.size() << " rows, ";
    if (!arr.empty()) {
        std::cout << arr[0].size() << " columns." << std::endl;
    } else {
        std::cout << "0 columns." << std::endl;
    }
}
std::vector<int> n_head_to_array(const std::string& input_str);
std::vector<int> f_head_to_array(const std::string& input_str);
std::vector<std::vector<Scalar>> nodes_coordinate_to_array(const std::string& input_str, int dim);
std::vector<std::vector<int>> face_thread_to_array(const std::string& input_str, int dim);
void print_2D_vector(const std::vector<std::vector<int>>& vec);
std::vector<std::vector<double>> readCoordinates(const std::string& filename);
std::vector<std::vector<int>> readFace(const std::string& filename);
// Mesh 类声明
class Mesh {
public:
    Mesh();

    void add_point(const std::vector<Scalar>& coordinates);

    void add_face(const std::vector<int>& face_data);

    // 查找包含指定点的所有面
    std::vector<std::vector<int>> find_faces_by_point_id(int point_id) const;

    // 解析 .msh 文件以提取网格数据
    void parse_msh_file(const std::string& file_path);

    // 获取点的坐标值
    const std::vector<std::vector<Scalar>>& get_points() const;

    // 获取面的信息
    const std::vector<std::vector<int>>& get_faces() const;

    

    // 重载 << 运算符，打印 Mesh 信息
    friend std::ostream& operator<<(std::ostream& os, const Mesh& mesh);

private:
    std::vector<std::vector<Scalar>> points; // 存储网格点坐标
    std::vector<std::vector<int>> faces;     // 存储面数据
   
};

// 点(三维向量)类定义
class Point {
public:
    Point();  // 默认构造函数，初始化为 (0, 0, 0)
    Point(Scalar x, Scalar y, Scalar z);  // 构造函数，接受三个坐标参数

    Scalar getX() const;
    Scalar getY() const;
    Scalar getZ() const;

    Point operator+(const Point& other) const;  // 向量加法
    Point operator-(const Point& other) const;  // 向量减法

    Scalar dot(const Point& other) const;  // 点乘
    Point cross(const Point& other) const;  // 叉乘

    Point operator*(Scalar scalar) const;  // 向量数乘
    friend Point operator*(Scalar scalar, const Point& point);  // 向量数乘（反向）

    Point operator/(Scalar scalar) const;  // 向量除以标量

    Point abs() const;  // 标准化

    Scalar magnitude() const;  // 计算向量的模（到原点的距离）
    Scalar distance(const Point& other) const;  // 计算两点之间的距离

    void print() const;  // 打印坐标

 

private:
    std::vector<Scalar> coordinates; // 存储坐标的向量
};

// Nodes 类，包含多个 Point 对象
class Nodes {
public:
    Nodes(const std::vector<std::vector<Scalar>>& coords);  // 构造函数，接收一个 n*3 的 vector<vector<Scalar>>

    size_t size() const;  // 获取节点数量
    const Point& getPoint(size_t index) const;  // 获取某个节点
    void print() const;  // 打印所有点
    const Point& operator[](size_t index) const;
private:
    std::vector<Point> points; // 存储所有节点的容器
};

// Face 类，定义面
class Face {
public:
    Face(const std::vector<int>& face_data, const Nodes& nodes);  // 构造函数，接收包含面数据的 vector

    void calculateNormal(const Nodes& nodes);  // 计算外法向量
    Point getNormal() const;  // 获取外法向量
    void printNormal() const;  // 打印法向量

    void calculateCenter(const Nodes& nodes);  // 计算面中心点
    Point getCenter() const;  // 获取中心点
    void printCenter() const;  // 打印中心点

    void calculateArea(const Nodes& nodes);  // 计算面面积
    Scalar getArea() const;  // 获取面面积
    void printArea() const;  // 打印面面积

private:
    int n0, n1, n2, n3;  // 四个顶点的索引
    int c0, c1;  // owner 和 neighbor 单元
    int zoneid, bctype, facetype;  // zone, boundary type, face type
    std::vector<int> nodes_indices;  // 存储面上所有的点索引
    Point normal;  // 面的外法向量
    Point center;  // 面的中心点
    Scalar area;   // 面的面积
};
//Faces类
class Faces {
public:
    Faces(const std::vector<std::vector<int>>& faces_data, const Nodes& nodes);

    void calculateAllNormals(const Nodes& nodes);
    void calculateAllCenters(const Nodes& nodes);
    void calculateAllAreas(const Nodes& nodes);
    

private:
    std::vector<Face> faces;
};

#endif // CFDMATH_H
