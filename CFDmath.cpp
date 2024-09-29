#include "CFDmath.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <regex>
#include <cmath>


// n_head_to_array 实现
std::vector<int> n_head_to_array(const std::string& input_str) {
    std::stringstream ss(input_str);
    std::vector<int> result;
    std::string element;
    int value;
    while (ss >> element) {
        if (result.size() == 2 || result.size() == 4) {
            value = std::stoi(element, nullptr, 16);  // 第3个元素以16进制解析
        } else {
            value = std::stoi(element);
        }
        result.push_back(value);
    }
    if (result.size() < 5) {
        throw std::runtime_error("cells/nodes 数组的元素数应至少为6");
    }
    return result;
}

// f_head_to_array 实现
std::vector<int> f_head_to_array(const std::string& input_str) {
    std::stringstream ss(input_str);
    std::vector<int> result;
    std::string element;
    while (ss >> element) {
        result.push_back(std::stoi(element, nullptr, 16));  // 全部为16进制解析
    }
    if (result.size() < 5) {
        throw std::runtime_error("cells/nodes 数组的元素数应至少为6");
    }
    return result;
}

// nodes_coordinate_to_array 实现
std::vector<std::vector<Scalar>> nodes_coordinate_to_array(const std::string& input_str, int dim) {
    std::vector<std::vector<Scalar>> result;
    std::stringstream ss(input_str);
    std::vector<Scalar> temp;
    Scalar value;

    while (ss >> value) {
        temp.push_back(value);
        if (temp.size() == dim) {
            result.push_back(temp);
            temp.clear();
        }
    }
    if (!temp.empty()) {
        throw std::runtime_error("元素坐标不匹配维度");
    }
    return result;
}

// face_thread_to_array 实现
std::vector<std::vector<int>> face_thread_to_array(const std::string& input_str, int dim) {
    std::vector<std::vector<int>> result;
    std::stringstream ss(input_str);
    std::vector<int> temp;
    int value;

    while (ss >> std::hex >> value) {
        temp.push_back(value);
        if (temp.size() == dim) {
            result.push_back(temp);
            temp.clear();
        }
    }
    return result;
}

// print_2D_vector 实现
void print_2D_vector(const std::vector<std::vector<int>>& vec) {
    for (const auto& row : vec) {
        for (const auto& element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}
// 函数：读取坐标，返回 n x 3 的坐标向量
std::vector<std::vector<double>> readCoordinates(const std::string& filename) {
    std::ifstream infile(filename);  // 打开输入文件
    if (!infile.is_open()) {  // 检查文件是否成功打开
        std::cerr << "打开文件出错: " << filename << std::endl;
        return {};  // 返回空向量
    }

    std::string line;  // 用于存储当前行内容
    bool readingCoordinates = false;  // 标记是否正在读取坐标
    std::vector<std::vector<double>> coordinates;  // 存储读取的坐标

    while (std::getline(infile, line)) {  // 逐行读取文件
        // 检查行是否匹配模式 "(10 ()("
        if (line.find("(10 (") != std::string::npos && line.find(")(") != std::string::npos) {
            readingCoordinates = true;  // 开始读取坐标
            continue;  // 继续到下一行
        }

        // 如果处于读取模式，检查坐标
        if (readingCoordinates) {
            // 检查当前行是否以 "))" 结束
            if (line.find("))") != std::string::npos) {
                break;  // 遇到结束，停止读取
            }

            // 从当前行读取三个浮点数
            std::istringstream iss(line);  // 创建字符串流
            double x, y, z;  // 定义坐标变量
            if (iss >> x >> y >> z) {  // 读取坐标
                coordinates.push_back({x, y, z});  // 存储坐标
            }
        }
    }

    return coordinates;  // 返回收集到的坐标
}


// 函数：读取十六进制面线程并转换为十进制整数，返回完整的面二维数组
std::vector<std::vector<int>> readFace(const std::string& filename) {
    std::ifstream infile(filename);  // 打开输入文件
    if (!infile.is_open()) {  // 检查文件是否成功打开
        std::cerr << "打开文件出错: " << filename << std::endl;
        return {};  // 返回空的二维数组
    }

    std::vector<std::vector<int>> allCombinedArrays;  // 存储所有结构的组合后的二维数组
    std::string line;  // 用于存储当前行内容

    while (std::getline(infile, line)) {  // 逐行读取文件
        // 检查行是否匹配模式 "(13 ("
        if (line.find("(13 (") != std::string::npos && line.find(")(") != std::string::npos) {
            std::vector<int> b;  // 存储第一个括号内的十进制数
            bool readingHex = true;  // 开始读取十六进制数

            std::size_t firstParenEnd = line.find(')');  // 查找第一个括号结束位置
            std::string firstBracketContent = line.substr(5, firstParenEnd - 5);  // 提取第一个括号内的内容

            // 转换第一个括号内的十六进制数
            std::istringstream iss(firstBracketContent);  // 创建字符串流
            std::string hexStr;  // 存储十六进制数的字符串
            while (iss >> hexStr) {  // 逐个读取十六进制数
                int decimalValue = std::stoi(hexStr, nullptr, 16);  // 转换为十进制整数
                b.push_back(decimalValue);  // 存储转换后的值
            }

            // 准备读取第二个括号内的内容
            std::vector<std::vector<int>> combinedArray;  // 存储当前结构的二维数组

            // 继续读取第二个括号内的内容
            while (std::getline(infile, line)) {  // 逐行读取文件
                // 检查当前行是否以 "))" 结束
                if (line.find("))") != std::string::npos) {
                    readingHex = false;  // 遇到结束，停止读取
                    break;  // 跳出循环
                }

                // 从当前行读取十六进制数
                std::istringstream iss(line);  // 创建字符串流
                std::vector<int> rowValues;  // 存储当前行的十进制值

                while (iss >> hexStr) {  // 逐个读取十六进制数
                    int decimalValue = std::stoi(hexStr, nullptr, 16);  // 转换为十进制整数
                    rowValues.push_back(decimalValue);  // 存储转换后的值
                }

                // 将当前行的值添加到二维数组中，并在前面加上数组 b
                if (!rowValues.empty()) {
                    rowValues.insert(rowValues.begin(), b.begin(), b.end());  // 在行开头插入数组 b
                    combinedArray.push_back(rowValues);  // 将行添加到当前结构的数组中
                }
            }

            // 将当前结构的二维数组添加到总数组中
            allCombinedArrays.insert(allCombinedArrays.end(), combinedArray.begin(), combinedArray.end());
        }
    }

    // 删除第二列、第三列和第五列
    for (auto& row : allCombinedArrays) {
        if (row.size() >= 5) {
            row.erase(row.begin() + 1);  // 删除第二列
            row.erase(row.begin() + 1);  // 删除第三列（现在的索引为1）
            if (row.size() >= 5) {
                row.erase(row.begin() + 3);  // 删除第五列（现在的索引为3）
            }
        }
    }

    return allCombinedArrays;  // 返回完整的二维数组
}
// Mesh 类的构造函数实现
Mesh::Mesh() :  points(), faces() {}


// 添加点
void Mesh::add_point(const std::vector<Scalar>& coordinates) {
    points.push_back(coordinates);
}

// 添加面
void Mesh::add_face(const std::vector<int>& face_data) {
    faces.push_back(face_data);
}

// 查找包含指定点的所有面
std::vector<std::vector<int>> Mesh::find_faces_by_point_id(int point_id) const {
    std::vector<std::vector<int>> faces_with_point;
    for (const auto& face : faces) {
        if (std::find(face.begin(), face.begin() + 4, point_id) != face.begin() + 4) {
            faces_with_point.push_back(face);
        }
    }
    return faces_with_point;
}


// 获取点的坐标值
const std::vector<std::vector<Scalar>>& Mesh::get_points() const {
    return points;
}


void Mesh::parse_msh_file(const std::string& file_path) {
    std::vector<std::vector<double>> coords = readCoordinates(file_path);  // 调用读取坐标的函数
    std::vector<std::vector<int>> facesThread = readFace(file_path);
 
    
    points =coords;
    faces = facesThread;
}


// 获取面的信息
const std::vector<std::vector<int>>& Mesh::get_faces() const {
    return faces;
}

// 重载 << 运算符，打印 Mesh 信息
std::ostream& operator<<(std::ostream& os, const Mesh& mesh) {
    os << "Points: " << std::endl;
    for (const auto& point : mesh.points) {
        for (auto coord : point) {
            os << coord << " ";
        }
        os << std::endl;
    }
    os << "Faces: " << std::endl;
    for (const auto& face : mesh.faces) {
        for (auto element : face) {
            os << element << " ";
        }
        os << std::endl;
    }
    return os;
}

// Point 类实现
Point::Point() : coordinates{0.0, 0.0, 0.0} {}
Point::Point(Scalar x, Scalar y, Scalar z) : coordinates{x, y, z} {}

Scalar Point::getX() const { return coordinates[0]; }
Scalar Point::getY() const { return coordinates[1]; }
Scalar Point::getZ() const { return coordinates[2]; }

Point Point::operator+(const Point& other) const {
    return Point(coordinates[0] + other.getX(),
                 coordinates[1] + other.getY(),
                 coordinates[2] + other.getZ());
}

Point Point::operator-(const Point& other) const {
    return Point(coordinates[0] - other.getX(),
                 coordinates[1] - other.getY(),
                 coordinates[2] - other.getZ());
}

Point Point::operator*(Scalar scalar) const {
    return Point(coordinates[0] * scalar,
                 coordinates[1] * scalar,
                 coordinates[2] * scalar);
}

Point Point::operator/(Scalar scalar) const {
    if (scalar == 0) {
        throw std::runtime_error("除以零错误");
    }
    return Point(coordinates[0] / scalar,
                 coordinates[1] / scalar,
                 coordinates[2] / scalar);
}

Scalar Point::dot(const Point& other) const {
    return coordinates[0] * other.getX() +
           coordinates[1] * other.getY() +
           coordinates[2] * other.getZ();
}

Point Point::cross(const Point& other) const {
    return Point(coordinates[1] * other.getZ() - coordinates[2] * other.getY(),
                 coordinates[2] * other.getX() - coordinates[0] * other.getZ(),
                 coordinates[0] * other.getY() - coordinates[1] * other.getX());
}

Scalar Point::magnitude() const {
    return std::sqrt(dot(*this));
}



// abs 实现
Point Point::abs() const {
    return Point(std::abs(coordinates[0]), std::abs(coordinates[1]), std::abs(coordinates[2]));
}


// 构造函数，接收一个 n*3 的 vector<vector<Scalar>>
Nodes::Nodes(const std::vector<std::vector<Scalar>>& coords) {
    for (const auto& coord : coords) {
        if (coord.size() != 3) {
            throw std::invalid_argument("每个坐标必须包含三个分量 (x, y, z)。");
        }
        // 创建 Point 对象并添加到 points 容器中
        points.emplace_back(coord[0], coord[1], coord[2]);
    }
}

// 获取节点数量
size_t Nodes::size() const {
    return points.size();
}

// 获取某个节点
const Point& Nodes::getPoint(size_t index) const {
  
    return points[index];
}

// 打印所有点
void Nodes::print() const {
    for (size_t i = 0; i < points.size(); ++i) {
        const Point& p = points[i];
        std::cout << "点 " << i << ": (" << p.getX() << ", " << p.getY() << ", " << p.getZ() << ")" << std::endl;
    }
}

// 重载运算符 [] 以访问分量
const Point& Nodes::operator[](size_t index) const {
    return getPoint(index); // 使用 getPoint 方法以安全方式访问点
}

// 构造函数，接收包含面数据的 vector
Face::Face(const std::vector<int>& face_data, const Nodes& nodes)
    : n0(face_data[3]), n1(face_data[4]), n2(face_data[5]), n3(face_data[6]),
      c0(face_data[7]), c1(face_data[8]),
      zoneid(face_data[0]), bctype(face_data[1]), facetype(face_data[2]) {
    // 存储面上所有的点索引
    nodes_indices = {n0, n1, n2, n3};
    calculateNormal(nodes);  // 计算外法向量
    calculateCenter(nodes);   // 计算面中心点
}

// 计算外法向量
void Face::calculateNormal(const Nodes& nodes) {
    // 通过右手法则计算法向量
    Point v1 = nodes[n1] - nodes[n0]; // n0 到 n1 的向量
    Point v2 = nodes[n2] - nodes[n0]; // n0 到 n2 的向量

    normal = v1.cross(v2);  // 计算法向量
    normal = normal.abs();   // 归一化
}

// 获取外法向量
Point Face::getNormal() const {
    return normal;
}

// 打印法向量
void Face::printNormal() const {
    Point n = getNormal();
    std::cout << "法向量: (" << n.getX() << ", " << n.getY() << ", " << n.getZ() << ")" << std::endl;
}

// 计算面中心点
void Face::calculateCenter(const Nodes& nodes) {
    // 计算面上四个顶点的平均值
    Point p0 = nodes[n0];
    Point p1 = nodes[n1];
    Point p2 = nodes[n2];
    Point p3 = nodes[n3];

    center = (p0 + p1 + p2 + p3) * (1.0 / 4.0);  // 取平均
}

// 获取中心点
Point Face::getCenter() const {
    return center;
}

// 打印中心点
void Face::printCenter() const {
    Point c = getCenter();
    std::cout << "中心点: (" << c.getX() << ", " << c.getY() << ", " << c.getZ() << ")" << std::endl;
}

// 计算面面积的方法
void Face::calculateArea(const Nodes& nodes) {
    // 获取四个顶点的坐标
    Point p0 = nodes.getPoint(n0);
    Point p1 = nodes.getPoint(n1);
    Point p2 = nodes.getPoint(n2);
    Point p3 = nodes.getPoint(n3);

    // 使用三角形的面积公式计算面面积
    double area1 = 0.5 * ((p1 - p0).cross(p2 - p0)).magnitude();
    double area2 = 0.5 * ((p2 - p0).cross(p3 - p0)).magnitude();
    
    area = area1 + area2;  // 总面积为两个三角形的面积之和
    
}

// 获取面面积的方法
Scalar Face::getArea() const {
    return area;
}

// 打印面面积的方法
void Face::printArea() const {
    std::cout << "Face Area: " << area << std::endl;
}
//构造函数 初始化Faces
Faces::Faces(const std::vector<std::vector<int>>& faces_data, const Nodes& nodes) {
    for (const auto& face_data : faces_data) {
        faces.emplace_back(face_data, nodes);
    }
}
//计算面法向量
void Faces::calculateAllNormals(const Nodes& nodes) {
    for (auto& face : faces) {
        face.calculateNormal(nodes);
    }
}
//计算面中心
void Faces::calculateAllCenters(const Nodes& nodes) {
    for (auto& face : faces) {
        face.calculateCenter(nodes);
    }
}
//计算所有面面积
void Faces::calculateAllAreas(const Nodes& nodes) {
    for (auto& face : faces) {
        face.calculateCenter(nodes);
    }
}

