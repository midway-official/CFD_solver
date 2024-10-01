#include "CFDmath.h"

void parseNodes(const std::string &filename, std::vector<std::vector<double>> &coordinates) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("无法打开文件");
    }

    std::string line;
    bool readingCoordinates = false; // 是否正在读取坐标

    while (std::getline(file, line)) {
        // 去除行首和行尾的空白字符
        line.erase(0, line.find_first_not_of(" \n\r\t")); // 左侧去空白
        line.erase(line.find_last_not_of(" \n\r\t") + 1); // 右侧去空白

        // 检查行是否以 '(' 结尾且以 '(10' 开头
        if (line.back() == '(') {
            if (line.find("(10") == 0) {
                readingCoordinates = true; // 开始读取坐标
                continue; // 跳过当前行
            }
        }

        // 如果正在读取坐标
        if (readingCoordinates) {
            // 检查行是否以 '))' 结尾
            if (line.find("))") != std::string::npos) {
                // 处理行中最后的坐标
                std::istringstream iss(line);
                double x, y, z;

                // 读取所有坐标，直到找到 '))'
                while (iss >> x >> y >> z) {
                    if (iss.eof()) break; // 如果到达行末，跳出
                    coordinates.push_back({x, y, z}); // 存储坐标
                }

                readingCoordinates = false; // 停止读取坐标
                continue; // 跳过当前行
            }

            // 从行中读取坐标
            std::istringstream iss(line);
            double x, y, z;

            while (iss >> x >> y >> z) {
                coordinates.push_back({x, y, z}); // 存储坐标
            }
        }
    }
}

void parseFaces(const std::string &filename, std::vector<std::vector<int>> &faces) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("无法打开文件");
    }

    std::string line;
    std::vector<int> headerInfo;
    bool readingFaces = false; // 是否正在读取面信息

    while (std::getline(file, line)) {
        // 去除行首和行尾的空白字符
        line.erase(0, line.find_first_not_of(" \n\r\t"));
        line.erase(line.find_last_not_of(" \n\r\t") + 1);

        // 检查行是否以 '(' 结尾
        if (line.back() == '(') {
            if (line.find("(13") == 0) {
                // 使用正则表达式提取头部信息
                std::regex headerRegex(R"(\(13\s+\((.*?)\))");
                std::smatch match;
                if (std::regex_search(line, match, headerRegex) && match.size() > 1) {
                    headerInfo.clear(); // 清空之前的头部信息
                    std::istringstream headerStream(match[1].str());
                    std::string hexValue;

                    while (headerStream >> hexValue) {
                        int decimalValue;
                        std::stringstream ss;
                        ss << std::hex << hexValue; // 将十六进制转换为十进制
                        ss >> decimalValue;
                        headerInfo.push_back(decimalValue); // 存储十进制值
                    }
                    readingFaces = true; // 开始读取面信息
                    continue; // 跳过当前行
                }
            }
        }

        // 如果正在读取面信息
        if (readingFaces) {
            // 检查行是否以 '))' 结尾
            if (line.find("))") != std::string::npos) {
                // 在停止之前处理当前行的面信息
                std::istringstream iss(line);
                std::string hexValue;
                std::vector<int> faceInfo;

                while (iss >> hexValue) {
                    int decimalValue;
                    std::stringstream ss;
                    ss << std::hex << hexValue; // 将十六进制转换为十进制
                    ss >> decimalValue;
                    faceInfo.push_back(decimalValue); // 存储面信息
                }

                // 合并头部信息和面信息
                faceInfo.insert(faceInfo.begin(), headerInfo.begin(), headerInfo.end());

                // 删除第 2、3 和 5 列
                if (faceInfo.size() >= 5) { // 确保面信息足够长
                    faceInfo.erase(faceInfo.begin() + 4); // 删除第 5 列
                    faceInfo.erase(faceInfo.begin() + 2); // 删除第 3 列
                    faceInfo.erase(faceInfo.begin() + 1); // 删除第 2 列
                }

                faces.push_back(faceInfo); // 添加面信息

                readingFaces = false; // 停止读取面信息
                continue; // 跳过当前行
            }

            // 从行中读取面信息
            std::istringstream iss(line);
            std::string hexValue;
            std::vector<int> faceInfo;

            while (iss >> hexValue) {
                int decimalValue;
                std::stringstream ss;
                ss << std::hex << hexValue; // 将十六进制转换为十进制
                ss >> decimalValue;
                faceInfo.push_back(decimalValue); // 存储面信息
            }

            // 合并头部信息和面信息
            faceInfo.insert(faceInfo.begin(), headerInfo.begin(), headerInfo.end());

            // 删除第 2、3 和 5 列
            if (faceInfo.size() >= 5) { // 确保面信息足够长
                faceInfo.erase(faceInfo.begin() + 4); // 删除第 5 列
                faceInfo.erase(faceInfo.begin() + 2); // 删除第 3 列
                faceInfo.erase(faceInfo.begin() + 1); // 删除第 2 列
            }

            faces.push_back(faceInfo); // 添加面信息
        }
    }
}


// 构造函数
Point::Point(const std::vector<double>& coordinates) {
    if (coordinates.size() != 3) {
        throw std::invalid_argument("Point must be initialized with a vector of size 3.");
    }
    this->coordinates = coordinates;
}

// 加法
Point Point::operator+(const Point& other) const {
    return Point({coordinates[0] + other.coordinates[0],
                   coordinates[1] + other.coordinates[1],
                   coordinates[2] + other.coordinates[2]});
}

// 减法
Point Point::operator-(const Point& other) const {
    return Point({coordinates[0] - other.coordinates[0],
                   coordinates[1] - other.coordinates[1],
                   coordinates[2] - other.coordinates[2]});
}

// 与double的乘法
Point Point::operator*(double scalar) const {
    return Point({coordinates[0] * scalar,
                   coordinates[1] * scalar,
                   coordinates[2] * scalar});
}

// 与double的除法
Point Point::operator/(double scalar) const {
    if (scalar == 0) {
        throw std::invalid_argument("Division by zero.");
    }
    return Point({coordinates[0] / scalar,
                   coordinates[1] / scalar,
                   coordinates[2] / scalar});
}

// 叉乘
Point Point::operator^(const Point& other) const {
    return Point({
        coordinates[1] * other.coordinates[2] - coordinates[2] * other.coordinates[1],
        coordinates[2] * other.coordinates[0] - coordinates[0] * other.coordinates[2],
        coordinates[0] * other.coordinates[1] - coordinates[1] * other.coordinates[0]
    });
}
//另一种叉乘
Point Point::cross(const Point& other) const{
     return Point({
        coordinates[1] * other.coordinates[2] - coordinates[2] * other.coordinates[1],
        coordinates[2] * other.coordinates[0] - coordinates[0] * other.coordinates[2],
        coordinates[0] * other.coordinates[1] - coordinates[1] * other.coordinates[0]
    });
}
//点乘
Scalar Point::dot(const Point& other) const{
     return 
        coordinates[1] * other.coordinates[1] + coordinates[2] * other.coordinates[2]+
        coordinates[0] * other.coordinates[0];
   
}
// 计算模长
double Point::magnitude() const {
    return std::sqrt(coordinates[0] * coordinates[0] +
                     coordinates[1] * coordinates[1] +
                     coordinates[2] * coordinates[2]);
}

// 单位化向量
Point Point::normalize() const {
    double mag = magnitude();
    if (mag == 0) {
        throw std::invalid_argument("Cannot normalize a zero vector.");
    }
    return *this / mag;
}


// 打印点坐标
void Point::print() const {
    std::cout << "(" << coordinates[0] << ", " << coordinates[1] << ", " << coordinates[2] << ")\n";
}


// 矢量场构造函数，接受一个 Point 向量，初始化为点（矢量）场
Field::Field(const std::vector<Point>& points) : isPointField(true), points(points) {}

//标量场 构造函数，接受一个标量向量，初始化为标量场
Field::Field(const std::vector<Scalar>& scalars) : isPointField(false), scalars(scalars) {}

// 重载加法运算符，返回两个 Field 对象的和
Field Field::operator+(const Field& other) const {
    checkCompatibility(other);
    if (isPointField) {//矢量场加法
        std::vector<Point> result;
        for (size_t i = 0; i < points.size(); ++i) {
            result.push_back(points[i] + other.points[i]);
        }
        return Field(result);
    } else {//标量场加法
        std::vector<Scalar> result;
        for (size_t i = 0; i < scalars.size(); ++i) {
            result.push_back(scalars[i] + other.scalars[i]);
        }
        return Field(result);
    }
}

// 重载减法运算符，返回两个 Field 对象的差
Field Field::operator-(const Field& other) const {
    checkCompatibility(other);
    if (isPointField) {//矢量场减法
        std::vector<Point> result;
        for (size_t i = 0; i < points.size(); ++i) {
            result.push_back(points[i] - other.points[i]);
        }
        return Field(result);
    } else {//标量场减法
        std::vector<Scalar> result;
        for (size_t i = 0; i < scalars.size(); ++i) {
            result.push_back(scalars[i] - other.scalars[i]);
        }
        return Field(result);
    }
}

// 重载乘法运算符，返回当前 Field 对象与标量的乘积
Field Field::operator*(double scalar) const {
    if (isPointField) {
        std::vector<Point> result;
        for (const auto& point : points) {
            result.push_back(point * scalar);
        }
        return Field(result);
    } else {
        std::vector<Scalar> result;
        for (const auto& value : scalars) {
            result.push_back(value * scalar);
        }
        return Field(result);
    }
}


// 重载乘法运算符，允许标量场与矢量场相乘
Field Field::operator*(const Field& other) const {
    if (isPointField && !other.isPointField) {
        // 当前是点场，另一个是标量场
        if (other.scalars.size() != points.size()) {
            throw std::invalid_argument("Size mismatch between point field and scalar field.");
        }
        std::vector<Point> result;
        for (size_t i = 0; i < points.size(); ++i) {
            result.push_back(points[i] * other.scalars[i]); // 矢量场与标量场相乘
        }
        return Field(result);
    } else if (!isPointField && other.isPointField) {
        // 当前是标量场，另一个是点场
        if (scalars.size() != other.points.size()) {
            throw std::invalid_argument("Size mismatch between scalar field and point field.");
        }
        std::vector<Point> result;
        for (size_t i = 0; i < other.points.size(); ++i) {
            result.push_back(other.points[i] * scalars[i]); // 矢量场与标量场相乘
        }
        return Field(result);
    } else {
        throw std::invalid_argument("Invalid field multiplication. One must be a scalar field and the other must be a point field.");
    }
}

// 矢量场之间的点乘
Field Field::dot(const Field& other) const {
    checkCompatibility(other);
    if (!isPointField || !other.isPointField) {
        throw std::invalid_argument("Dot product can only be computed between point fields.");
    }
    std::vector<Scalar> result;
    for (size_t i = 0; i < points.size(); ++i) {
        result.push_back(points[i].dot(other.points[i])); // 计算点乘
    }
    return Field(result);
}

// 矢量场之间的叉乘
Field Field::cross(const Field& other) const {
    checkCompatibility(other);
    if (!isPointField || !other.isPointField) {
        throw std::invalid_argument("Cross product can only be computed between point fields.");
    }
    std::vector<Point> result;
    for (size_t i = 0; i < points.size(); ++i) {
        result.push_back(points[i] ^ other.points[i]); // 计算叉乘
    }
    return Field(result);
}

// 计算当前点场的模长，返回标量场
Field Field::magnitude() const {
    if (!isPointField) {
        throw std::invalid_argument("Magnitude can only be computed for point fields.");
    }
    std::vector<Scalar> magnitudes;
    for (const auto& point : points) {
        magnitudes.push_back(point.magnitude());
    }
    return Field(magnitudes);
}

// 标准化当前的点场，返回标准化后的点场
Field Field::normalize() const {
    if (!isPointField) {
        throw std::invalid_argument("Normalization can only be performed on point fields.");
    }
    std::vector<Point> normalizedPoints;
    for (const auto& point : points) {
        normalizedPoints.push_back(point.normalize());
    }
    return Field(normalizedPoints);
}

// 打印当前 Field 的内容
void Field::print() const {
    if (isPointField) {
        for (const auto& point : points) {
            point.print();
        }
    } else {
        for (const auto& scalar : scalars) {
            std::cout << scalar << " ";
        }
        std::cout << std::endl;
    }
}

// 检查两个 Field 是否兼容
void Field::checkCompatibility(const Field& other) const {
    if (isPointField != other.isPointField || 
        (isPointField && points.size() != other.points.size()) ||
        (!isPointField && scalars.size() != other.scalars.size())) {
        throw std::invalid_argument("Fields are not compatible for operation.");
    }
}

// 面集构造函数
Faces::Faces(const std::vector<std::vector<double>>& nodeCoordinates,
             const std::vector<std::vector<int>>& faceInfo) 
    : nodes(nodeCoordinates), faces(faceInfo) {
    // 检查节点坐标的有效性
    for (const auto& node : nodes) {
        if (node.size() != 3) {
            throw std::invalid_argument("Each node must have three coordinates.");
        }
    }

    // 检查面的信息的有效性
    for (const auto& face : faces) {
        if (face.size() != 8) {
            throw std::invalid_argument("Each face info must have eight elements (zoneid, bctype, n0, n1, n2, n3, c0, c1).");
        }
        // 检查节点编号的有效性（从1开始）
        for (int i = 2; i < 6; ++i) { // n0 到 n3
            if (face[i] < 1 || face[i] -1 > nodes.size()) {
                throw std::invalid_argument("Node indices must be within the range of available nodes.");
            }
        }
    }
}



// 打印面的信息
void Faces::print() const {
    for (const auto& face : faces) {
        std::cout << "ZoneID: " << face[0]
                  << ", BCType: " << face[1]
                  << ", Nodes: (" << face[2] << ", " << face[3] << ", "
                  << face[4] << ", " << face[5] << ")"
                  << ", C0: " << face[6]
                  << ", C1: " << face[7] << "\n";
    }
}
//计算面中心点，返回一个矢量场
Field Faces::calculateAllCenters() const {
    std::vector<Point> centers;

    for (const auto& face : faces) {
        // 确保面中包含足够的节点
        if (face.size() < 6) {
            throw std::runtime_error("Face does not have enough nodes.");
        }

        double x = 0.0, y = 0.0, z = 0.0;

        for (int i = 2; i <= 5; ++i) {
            int nodeIndex = face[i] -1 ; // 调整索引
            if (nodeIndex < 0 || nodeIndex >= nodes.size()) {
                std::cerr << "Error: Node index " << nodeIndex << " out of bounds. "
                          << "nodes.size() = " << nodes.size() << std::endl;
                throw std::runtime_error("Node index out of bounds.");
            }

            x += nodes[nodeIndex][0];
            y += nodes[nodeIndex][1];
            z += nodes[nodeIndex][2];
        }

        // 计算中心点并添加到结果中
        centers.emplace_back(std::vector<double>{x / 4.0, y / 4.0, z / 4.0});
    }

    return Field(centers);
}

//计算面法向量 由c1指向c0 n0 n1 n2确定
Field Faces::calculateAllNormals() const {
    std::vector<Point> normals;

    for (const auto& face : faces) {
        // 确保面中包含足够的节点
        

        // 获取节点坐标
        int nodeIndex0 = face[2] - 1; // n0
        int nodeIndex1 = face[3] - 1; // n1
        int nodeIndex2 = face[4] - 1; // n2

        if (nodeIndex0 < 0 || nodeIndex0 >= nodes.size() ||
            nodeIndex1 < 0 || nodeIndex1 >= nodes.size() ||
            nodeIndex2 < 0 || nodeIndex2 >= nodes.size()) {
            std::cerr << "Error: Node index out of bounds." << std::endl;
            throw std::runtime_error("Node index out of bounds.");
        }

        // 计算边向量
        std::vector<double> v1 = {
            nodes[nodeIndex1][0] - nodes[nodeIndex0][0],
            nodes[nodeIndex1][1] - nodes[nodeIndex0][1],
            nodes[nodeIndex1][2] - nodes[nodeIndex0][2]
        };

        std::vector<double> v2 = {
            nodes[nodeIndex2][0] - nodes[nodeIndex0][0],
            nodes[nodeIndex2][1] - nodes[nodeIndex0][1],
            nodes[nodeIndex2][2] - nodes[nodeIndex0][2]
        };

        // 计算法向量（叉乘）
        double nx = v1[1] * v2[2] - v1[2] * v2[1];
        double ny = v1[2] * v2[0] - v1[0] * v2[2];
        double nz = v1[0] * v2[1] - v1[1] * v2[0];

        // 法向量归一化
        double length = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (length > 0) {
            nx /= length;
            ny /= length;
            nz /= length;
        }

        // 将法向量添加到结果中
        normals.emplace_back(std::vector<double>{nx, ny, nz});
    }

    return Field(normals);
}

Field Faces::calculateAllAreas() const {
    std::vector<Scalar> areas;

    for (const auto& face : faces) {
        // 确保面中包含足够的节点
        if (face.size() < 4) { // 至少需要四个节点
            throw std::runtime_error("Face does not have enough nodes.");
        }

        // 获取节点坐标
        int nodeIndex0 = face[2] - 1; // n0
        int nodeIndex1 = face[3] - 1; // n1
        int nodeIndex2 = face[4] - 1; // n2
        int nodeIndex3 = face[5] - 1; // n3

        if (nodeIndex0 < 0 || nodeIndex0 >= nodes.size() ||
            nodeIndex1 < 0 || nodeIndex1 >= nodes.size() ||
            nodeIndex2 < 0 || nodeIndex2 >= nodes.size() ||
            nodeIndex3 < 0 || nodeIndex3 >= nodes.size()) {
            std::cerr << "Error: Node index out of bounds." << std::endl;
            throw std::runtime_error("Node index out of bounds.");
        }

        // 计算边向量
        std::vector<double> v1 = {
            nodes[nodeIndex1][0] - nodes[nodeIndex0][0],
            nodes[nodeIndex1][1] - nodes[nodeIndex0][1],
            nodes[nodeIndex1][2] - nodes[nodeIndex0][2]
        };

        std::vector<double> v2 = {
            nodes[nodeIndex2][0] - nodes[nodeIndex0][0],
            nodes[nodeIndex2][1] - nodes[nodeIndex0][1],
            nodes[nodeIndex2][2] - nodes[nodeIndex0][2]
        };

        std::vector<double> v3 = {
            nodes[nodeIndex3][0] - nodes[nodeIndex0][0],
            nodes[nodeIndex3][1] - nodes[nodeIndex0][1],
            nodes[nodeIndex3][2] - nodes[nodeIndex0][2]
        };

        // 计算法向量（叉乘）并计算两个三角形的面积
        double nx1 = v1[1] * v2[2] - v1[2] * v2[1];
        double ny1 = v1[2] * v2[0] - v1[0] * v2[2];
        double nz1 = v1[0] * v2[1] - v1[1] * v2[0];
        double area1 = 0.5 * std::sqrt(nx1 * nx1 + ny1 * ny1 + nz1 * nz1);

        double nx2 = v2[1] * v3[2] - v2[2] * v3[1];
        double ny2 = v2[2] * v3[0] - v2[0] * v3[2];
        double nz2 = v2[0] * v3[1] - v2[1] * v3[0];
        double area2 = 0.5 * std::sqrt(nx2 * nx2 + ny2 * ny2 + nz2 * nz2);

        // 总面积
        double totalArea = area1 + area2;

        // 将面积添加到结果中
        areas.push_back(totalArea);
    }

    return Field(areas);
}