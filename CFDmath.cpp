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


// 构造函数
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
void Faces::printFaces() const {
    for (const auto& face : faces) {
        std::cout << "ZoneID: " << face[0]
                  << ", BCType: " << face[1]
                  << ", Nodes: (" << face[2] << ", " << face[3] << ", "
                  << face[4] << ", " << face[5] << ")"
                  << ", C0: " << face[6]
                  << ", C1: " << face[7] << "\n";
    }
}

std::vector<Point> Faces::calculateAllCenters() const {
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

    return centers;
}