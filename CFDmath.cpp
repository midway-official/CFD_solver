#include "CFDmath.h"

void parseNodes(const std::string& filename, std::vector<std::vector<double>>& coordinates) {
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
                    coordinates.push_back({ x, y, z }); // 存储坐标
                }

                readingCoordinates = false; // 停止读取坐标
                continue; // 跳过当前行
            }

            // 从行中读取坐标
            std::istringstream iss(line);
            double x, y, z;

            while (iss >> x >> y >> z) {
                coordinates.push_back({ x, y, z }); // 存储坐标
            }
        }
    }
}

void parseFaces(const std::string& filename, std::vector<std::vector<int>>& faces) {
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

void parseFacetree(const std::string& filename, std::vector<std::vector<int>>& facetree) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("无法打开文件");
    }

    std::string line;
    std::vector<int> headerInfo;
    bool readingFacetree = false; // 是否正在读取面tree信息

    while (std::getline(file, line)) {
        // 去除行首和行尾的空白字符
        line.erase(0, line.find_first_not_of(" \n\r\t"));
        line.erase(line.find_last_not_of(" \n\r\t") + 1);

        // 检查行是否以 '(' 结尾
        if (line.back() == '(') {
            if (line.find("(59") == 0) {
                std::regex headerRegex(R"(\(59\s+\((.*?)\))");
                std::smatch match;
                if (std::regex_search(line, match, headerRegex) && match.size() > 1) {
                    headerInfo.clear();
                    std::istringstream headerStream(match[1].str());
                    std::string hexValue;

                    while (headerStream >> hexValue) {
                        int decimalValue;
                        std::stringstream ss;
                        ss << std::hex << hexValue;
                        ss >> decimalValue;
                        headerInfo.push_back(decimalValue);
                    }
                    readingFacetree = true;
                    continue;
                }
            }
        }

        // 如果正在读取面tree信息
        if (readingFacetree) {
            if (line.find("))") != std::string::npos) {
                std::istringstream iss(line);
                std::string hexValue;
                std::vector<int> faceInfo;

                while (iss >> hexValue) {
                    int decimalValue;
                    std::stringstream ss;
                    ss << std::hex << hexValue;
                    ss >> decimalValue;
                    faceInfo.push_back(decimalValue);
                }

                // 只取最后两列
                std::vector<int> newRow;
                newRow.push_back(facetree.size() + 1); // 插入编号
                newRow.push_back(faceInfo[faceInfo.size() - 2]); // 倒数第二列
                newRow.push_back(faceInfo.back()); // 最后一列

                facetree.push_back(newRow);

                readingFacetree = false;
                continue;
            }

            std::istringstream iss(line);
            std::string hexValue;
            std::vector<int> faceInfo;

            while (iss >> hexValue) {
                int decimalValue;
                std::stringstream ss;
                ss << std::hex << hexValue;
                ss >> decimalValue;
                faceInfo.push_back(decimalValue);
            }

            // 只取最后两列
            std::vector<int> newRow;
            newRow.push_back(facetree.size() + 1); // 插入编号
            newRow.push_back(faceInfo[faceInfo.size() - 2]); // 倒数第二列
            newRow.push_back(faceInfo.back()); // 最后一列

            facetree.push_back(newRow);
        }
    }
}
//计算cell信息


std::vector<std::vector<int>> processCelldata(std::vector<std::vector<int>> faceData, const std::vector<std::vector<int>>& facetree) {
    int maxCells = 0;

    // 查找faceData中最后两列的最大值作为总cell数量
    for (const auto& row : faceData) {
        if (row.size() >= 2) {
            maxCells = std::max(maxCells, std::max(row[row.size() - 2], row[row.size() - 1]));
        }
    }

    std::vector<bool> toRemove(faceData.size(), false);
    std::vector<int> freeindex;

    // 读取面树信息的第二列和第三列
    for (const auto& row : facetree) {
        if (row.size() >= 3) {
            freeindex.push_back(row[1] - 1); // F1, 转换为索引
            freeindex.push_back(row[2] - 1); // F2, 转换为索引
        }
    }

    // 标记需要删除的行
    for (int index : freeindex) {
        if (index >= 0 && index < toRemove.size()) {
            toRemove[index] = true;
        }
    }

    // 使用拷贝进行删除标记为true的行
    std::vector<std::vector<int>> modifiedFaceData = faceData; // 创建拷贝
    modifiedFaceData.erase(std::remove_if(modifiedFaceData.begin(), modifiedFaceData.end(),
        [&toRemove, &modifiedFaceData](const std::vector<int>& face) {
            int index = &face - &modifiedFaceData[0];  // 获取当前行的索引
            return index < modifiedFaceData.size() && toRemove[index]; // 判断是否需要删除
        }), modifiedFaceData.end());

    std::vector<std::vector<int>> cellData(maxCells);

    // 处理celldata
    for (int i = 1; i <= maxCells; ++i) {
        std::vector<int> index;
        std::unordered_set<int> resultSet;

        // 提取modifiedFaceData中最后两列有i的行的编号
        for (size_t row = 0; row < modifiedFaceData.size(); ++row) {
            if (modifiedFaceData[row][modifiedFaceData[row].size() - 2] == i || modifiedFaceData[row][modifiedFaceData[row].size() - 1] == i) {
                index.push_back(row + 1); // 从1开始的编号
            }
        }

        // 读取对应行并处理
        for (int idx : index) {
            const auto& row = modifiedFaceData[idx - 1]; // 转换为0开始的索引

            // 添加modifiedFaceData的行中的第3, 4, 5, 6个元素到结果集中
            for (size_t j = 2; j < 6 && j < row.size(); ++j) {
                resultSet.insert(row[j]);
            }
        }

        // 转换unordered_set为vector并存储到cellData
        cellData[i - 1] = std::vector<int>(resultSet.begin(), resultSet.end());
    }

    return cellData;
}

void printMaxValue(const std::vector<std::vector<int>>& vec) {
    if (vec.empty() || vec[0].empty()) {
        std::cout << "The vector is empty." << std::endl; // 如果向量为空
        return;
    }

    int maxValue = std::numeric_limits<int>::min(); // 初始化为最小整数

    // 遍历二维向量，找到最大值
    for (const auto& innerVec : vec) {
        for (int value : innerVec) {
            if (value > maxValue) {
                maxValue = value; // 更新最大值
            }
        }
    }

    std::cout << "Maximum value: " << maxValue << std::endl; // 打印最大值
}

void printNodeData(const std::vector<std::vector<double>>& nodeData) {
    std::cout << "Node Data:" << std::endl; // 打印节点数据标题
    for (const auto& node : nodeData) {
        for (double value : node) {
            std::cout << value << " "; // 打印每个节点的值
        }
        std::cout << std::endl; // 换行
    }
}

void printData(const std::vector<std::vector<int>>& faceData) {
    std::cout << "Face Data:" << std::endl; // 打印面数据标题
    size_t maxRows = std::min(faceData.size(), size_t(100000)); // 确保不超过100行
    for (size_t row = 0; row < maxRows; ++row) {
        for (int value : faceData[row]) {
            std::cout << value << " "; // 打印每个面的值
        }
        std::cout << std::endl; // 换行
    }
}

void printFaces(const std::vector<std::vector<int>>& faces) {
    // 打印每个面
    for (const auto& face : faces) {
        std::cout << "(";
        for (size_t i = 0; i < face.size(); ++i) {
            std::cout << face[i]; // 打印面中的每个值
            if (i < face.size() - 1) {
                std::cout << ", "; // 添加逗号分隔
            }
        }
        std::cout << ")" << std::endl; // 每个面换行
    }
}

void printMaxSubvectorLength(const std::vector<std::vector<int>>& vec) {
    int maxLength = 0; // 初始化最大长度
    // 遍历每个子向量，找到最大长度
    for (const auto& subvec : vec) {
        maxLength = std::max(maxLength, static_cast<int>(subvec.size())); // 更新最大长度
    }
    std::cout << maxLength << std::endl; // 打印最大子向量长度
}


// 向量加法函数定义
std::vector<double> addVectors(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    // 检查向量大小是否相同
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Error: Vectors must be of the same size.");
    }

    std::vector<double> result(vec1.size()); // 初始化结果向量

    // 按元素相加
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] + vec2[i];
    }

    return result; // 返回结果向量
}

// point构造函数

Point::Point(const std::vector<double>& coordinates) {
    if (coordinates.size() != 3) {
        throw std::invalid_argument("必须使用长度为3的矢量构造点");
    }
    this->coordinates = coordinates;
}

// 加法
Point Point::operator+(const Point& other) const {
    return Point({ coordinates[0] + other.coordinates[0],
                   coordinates[1] + other.coordinates[1],
                   coordinates[2] + other.coordinates[2] });
}

// 减法
Point Point::operator-(const Point& other) const {
    return Point({ coordinates[0] - other.coordinates[0],
                   coordinates[1] - other.coordinates[1],
                   coordinates[2] - other.coordinates[2] });
}

// 与double的乘法
Point Point::operator*(double scalar) const {
    return Point({ coordinates[0] * scalar,
                   coordinates[1] * scalar,
                   coordinates[2] * scalar });
}

// 与double的除法
Point Point::operator/(double scalar) const {
    if (scalar == 0) {
        throw std::invalid_argument("除以0.");
    }
    return Point({ coordinates[0] / scalar,
                   coordinates[1] / scalar,
                   coordinates[2] / scalar });
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
Point Point::cross(const Point& other) const {
    return Point({
       coordinates[1] * other.coordinates[2] - coordinates[2] * other.coordinates[1],
       coordinates[2] * other.coordinates[0] - coordinates[0] * other.coordinates[2],
       coordinates[0] * other.coordinates[1] - coordinates[1] * other.coordinates[0]
        });
}
//点乘
Scalar Point::dot(const Point& other) const {
    return
        coordinates[1] * other.coordinates[1] + coordinates[2] * other.coordinates[2] +
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
        throw std::invalid_argument("无法单位化零向量.");
    }
    return *this / mag;
}


// 打印点坐标
void Point::print() const {
    std::cout << "(" << coordinates[0] << ", " << coordinates[1] << ", " << coordinates[2] << ")\n";
}
// 获取坐标
double Point::getCoordinate(size_t index) const {
    if (index >= 3) throw std::out_of_range("Index out of range");
    return coordinates[index];
}

// 设置坐标
void Point::setCoordinate(size_t index, double value) {
    if (index >= 3) throw std::out_of_range("Index out of range");
    coordinates[index] = value;
}




// 矢量场构造函数，接受一个 Point 向量，初始化为点（矢量）场
Field::Field(const std::vector<Point>& points) : isPointField(true), points(points) {}

//标量场 构造函数，接受一个标量向量，初始化为标量场
Field::Field(const std::vector<Scalar>& scalars) : isPointField(false), scalars(scalars) {}
//长度和标量构造标量场
Field::Field(size_t elementCount, Scalar value) : isPointField(false) {
    scalars.resize(elementCount, value);
}
//长度和矢量构造矢量场
Field::Field(size_t elementCount, const std::vector<double>& coordinates) : isPointField(true) {
    for (size_t i = 0; i < elementCount; ++i) {
        points.emplace_back(coordinates);
    }
}
// 实现新增构造函数
Field::Field(const std::vector<std::vector<Scalar>>& vectorField) {
    if (vectorField.empty() || vectorField[0].size() != 3) {
        throw std::invalid_argument("Each element of the vectorField must have exactly three components.");
    }

    isPointField = true; // 设置为点场
    points.resize(vectorField.size());
    for (size_t i = 0; i < vectorField.size(); ++i) {
        points[i] = Point(vectorField[i]); // 初始化 Point 对象
    }
}
// 重载加法运算符，返回两个 Field 对象的和
Field Field::operator+(const Field& other) const {
    checkCompatibility(other);
    if (isPointField) {//矢量场加法
        std::vector<Point> result;
        for (size_t i = 0; i < points.size(); ++i) {
            result.push_back(points[i] + other.points[i]);
        }
        return Field(result);
    }
    else {//标量场加法
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
    }
    else {//标量场减法
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
    }
    else {
        std::vector<Scalar> result;
        for (const auto& value : scalars) {
            result.push_back(value * scalar);
        }
        return Field(result);
    }
}
// 标量与 Field 对象的乘法
Field operator*(double scalar, const Field& field) {
    return field * scalar; // 利用已经实现的 Field * scalar 的实现
}

// 重载乘法运算符，允许标量场与矢量场相乘
Field Field::operator*(const Field& other) const {
    if (isPointField && !other.isPointField) {
        // 当前是点场，另一个是标量场
        if (other.scalars.size() != points.size()) {
            throw std::invalid_argument("矢量场与标量场大小不相容.");
        }
        std::vector<Point> result;
        for (size_t i = 0; i < points.size(); ++i) {
            result.push_back(points[i] * other.scalars[i]); // 矢量场与标量场相乘
        }
        return Field(result);
    }
    else if (!isPointField && other.isPointField) {
        // 当前是标量场，另一个是点场
        if (scalars.size() != other.points.size()) {
            throw std::invalid_argument("矢量场与标量场大小不相容d.");
        }
        std::vector<Point> result;
        for (size_t i = 0; i < other.points.size(); ++i) {
            result.push_back(other.points[i] * scalars[i]); // 矢量场与标量场相乘
        }
        return Field(result);
    }
    else {
        throw std::invalid_argument("必须是一个矢量场一个标量场.");
    }
}

// 矢量场之间的点乘
Field Field::dot(const Field& other) const {
    checkCompatibility(other);
    if (!isPointField || !other.isPointField) {
        throw std::invalid_argument("只有矢量场之间可以点乘.");
    }
    std::vector<Scalar> result;
    for (size_t i = 0; i < points.size(); ++i) {
        result.push_back(points[i].dot(other.points[i])); // 计算点乘
    }
    return Field(result);
}
//场除标量
Field Field::operator/(double scalar) const {
    if (scalar == 0) {
        throw std::invalid_argument("Division by zero.");
    }

    // 创建结果对象，使用拷贝构造
    Field result(*this);

    if (isPointField) {
        // 如果是矢量场，逐个点进行除法
        for (Point& p : result.points) {
            p = p / scalar; // 调用 Point 的除法运算符
        }
    }
    else {
        // 如果是标量场，逐个标量进行除法
        for (Scalar& s : result.scalars) {
            s /= scalar; // 直接对标量进行除法
        }
    }

    return result;
}
Field Field::operator/(const Field& other) const {
    // 确保两个 Field 对象兼容
    if (size() != other.size()) {
        throw std::invalid_argument("两个 Field 对象大小不匹配");
    }

    // 根据被除的场类型初始化结果对象
    if (isPointField && !other.isPointField) {
        // 矢量场除以标量场，初始化结果为矢量场
        Field result(size(), std::vector<double>(3, 0.0)); // 初始化矢量场
        for (size_t i = 0; i < size(); ++i) {
            for (size_t j = 0; j < 3; ++j) { // 三维空间
                double newValue = points[i].getCoordinate(j) / other.scalarAt(i);
                result.pointAt(i).setCoordinate(j, newValue);
            }
        }
        return result; // 返回矢量场
    }
    else if (!isPointField && !other.isPointField) {
        // 标量场除以标量场，初始化结果为标量场
        Field result(size(), 0.0); // 初始化标量场
        for (size_t i = 0; i < size(); ++i) {
            result.scalarAt(i) = scalars[i] / other.scalarAt(i);
        }
        return result; // 返回标量场
    }
    else {
        throw std::invalid_argument("不允许标量场除以矢量场或两个矢量场进行除法");
    }
}


// 矢量场之间的叉乘
Field Field::cross(const Field& other) const {
    checkCompatibility(other);
    if (!isPointField || !other.isPointField) {
        throw std::invalid_argument("只有矢量场之间可以叉乘.");
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
        throw std::invalid_argument("只能返回矢量场的模长标量场.");
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
        throw std::invalid_argument("只能返回矢量场的标准化矢量场.");
    }
    std::vector<Point> normalizedPoints;
    for (const auto& point : points) {
        normalizedPoints.push_back(point.normalize());
    }
    return Field(normalizedPoints);
}

//索引访问场中某个点或标量从0开始
Point& Field::pointAt(size_t index) {
    return points[index];
}

const Point& Field::pointAt(size_t index) const {
    return points[index];
}

Scalar& Field::scalarAt(size_t index) {
    return scalars[index];
}

const Scalar& Field::scalarAt(size_t index) const {
    return scalars[index];
}
// 打印当前 Field 的内容
void Field::print() const {
    if (isPointField) {
        for (const auto& point : points) {
            point.print();
        }
    }
    else {
        for (const auto& scalar : scalars) {
            std::cout << scalar << " ";
        }
        std::cout << std::endl;
    }
}
//场元素个数
size_t Field::size() const {
    return isPointField ? points.size() : scalars.size();
}
// 获取标量场
std::vector<Scalar> Field::getScalarVector() const {
    return scalars; // 直接返回标量场数据
}

// 获取点场
std::vector<std::vector<Scalar>> Field::getPointVector() const {
    std::vector<std::vector<Scalar>> result;

    // 遍历 points，将每个点的坐标放入 result
    for (const auto& point : points) {
        result.push_back({point.getCoordinate(0), point.getCoordinate(1), point.getCoordinate(2)});
    }

    return result; // 返回点场数据
}
// 获取向量场的 x 分量
std::vector<double> Field::getXComponent() const {
    std::vector<double> xComponents;
    for (const auto& point : points) {
        xComponents.push_back(point.getCoordinate(0)); // 获取 x 坐标
    }
    return xComponents;
}

// 获取向量场的 y 分量
std::vector<double> Field::getYComponent() const {
    std::vector<double> yComponents;
    for (const auto& point : points) {
        yComponents.push_back(point.getCoordinate(1)); // 获取 y 坐标
    }
    return yComponents;
}

// 获取向量场的 z 分量
std::vector<double> Field::getZComponent() const {
    std::vector<double> zComponents;
    for (const auto& point : points) {
        zComponents.push_back(point.getCoordinate(2)); // 获取 z 坐标
    }
    return zComponents;
}

// 检查两个 Field 是否兼容
void Field::checkCompatibility(const Field& other) const {
    if (isPointField != other.isPointField ||
        (isPointField && points.size() != other.points.size()) ||
        (!isPointField && scalars.size() != other.scalars.size())) {
        throw std::invalid_argument("场类型不兼容，无法进行操作");
    }
}

// 面集构造函数
Faces::Faces() : nodes(), faces() {}
Faces::Faces(const std::vector<std::vector<double>>& nodeCoordinates,
    const std::vector<std::vector<int>>& faceInfo)
    : nodes(nodeCoordinates), faces(faceInfo) {
    // 检查节点坐标的有效性
    for (const auto& node : nodes) {
        if (node.size() != 3) {
            throw std::invalid_argument("每个点必须有3个坐标.");
        }
    }

    // 检查面的信息的有效性
    for (const auto& face : faces) {
        if (face.size() != 8) {
            throw std::invalid_argument("面必须包含8个信息 (zoneid, bctype, n0, n1, n2, n3, c0, c1).");
        }
        // 检查节点编号的有效性（从1开始）
        for (int i = 2; i < 6; ++i) { // n0 到 n3
            if (face[i] < 1 || face[i] - 1 > nodes.size()) {
                throw std::invalid_argument("节点必须在可用范围内.");
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
            int nodeIndex = face[i] - 1; // 调整索引
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
//计算所有面面积
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

// 访问索引为 i 的面实现
const std::vector<int>& Faces::getFace(size_t i) const {
    if (i >= faces.size()) {
        throw std::out_of_range("Face index out of range.");
    }
    return faces[i];
}

// 修改 faces 中第一个元素为 faceIndex 的行的第二个元素实现
void Faces::setBctypeForFace(int zoneIndex, int bctype) {
    for (auto& face : faces) {
        if (face[0] == zoneIndex) {
            face[1] = bctype; // 修改第二个元素
            return;
        }
    }
    throw std::runtime_error("Face with the specified index not found.");
}
//计算face个数
size_t Faces::size() const {
    return faces.size();
}
//修改场中的值
void Field::setScalar(size_t index, Scalar value) {
    if (!isPointField) {
        if (index < scalars.size()) {
            scalars[index] = value;
        }
        else {
            throw std::out_of_range("Index out of range for scalars.");
        }
    }
    else {
        std::cerr << "Warning: Attempting to modify a scalar in a point field.\n";
    }
}

void Field::setPoint(size_t index, const Point& point) {
    if (isPointField) {
        if (index < points.size()) {
            points[index] = point;
        }
        else {
            throw std::out_of_range("Index out of range for points.");
        }
    }
    else {
        std::cerr << "Warning: Attempting to modify a point in a scalar field.\n";
    }
}

// cell集构造函数
Cells::Cells() : nodes(), cells() {}
Cells::Cells(const std::vector<std::vector<double>>& nodeCoordinates,
    const std::vector<std::vector<int>>& cellsdata)
    : nodes(nodeCoordinates), cells(cellsdata) {
    // 检查节点坐标的有效性
    for (const auto& node : nodes) {
        if (node.size() != 3) {
            throw std::invalid_argument("每个点必须有3个坐标.");
        }
    }

    // 检查面的信息的有效性
    for (const auto& cell : cells) {
        if (cell.size() != 8) {
            throw std::invalid_argument("每个体必须包含8个电 .");

        }
        // 检查节点编号的有效性（从1开始）
        for (int i = 2; i < 6; ++i) { // n0 到 n3
            if (cell[i] < 1 || cell[i] - 1 > nodes.size()) {
                throw std::invalid_argument("节点必须在可用范围内.");
            }
        }
    }
}
//计算cell中心点，返回一个矢量场
Field Cells::calculateAllCenters() const {
    std::vector<Point> centers;

    for (const auto& cell : cells) {


        double x = 0.0, y = 0.0, z = 0.0;

        for (int i = 0; i <= 7; ++i) {
            int nodeIndex = cell[i] - 1; // 调整索引
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
        centers.emplace_back(std::vector<double>{x / 8.0, y / 8.0, z / 8.0});
    }

    return Field(centers);
}
//计算体积
Field Cells::calculateAllVolumes() const {
    std::vector<Scalar> volumes;

    for (const auto& cell : cells) {
        if (cell.size() < 8) {
            throw std::runtime_error("Cell does not have enough nodes.");
        }

        // 获取节点坐标
        std::vector<std::vector<double>> vertices(8);
        for (int i = 0; i < 8; ++i) {
            int nodeIndex = cell[i] - 1;
            if (nodeIndex < 0 || nodeIndex >= nodes.size()) {
                throw std::runtime_error("Node index out of bounds.");
            }
            vertices[i] = nodes[nodeIndex];
        }

        // 按z坐标排序
        std::sort(vertices.begin(), vertices.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
            return a[2] < b[2];
            });

        // 获取z最大值的后四个点
        std::vector<std::vector<double>> topVertices = { vertices[4], vertices[5], vertices[6], vertices[7] };

        // 按逆时针排序点
        std::vector<double> centroid = { 0.0, 0.0 };
        for (const auto& v : topVertices) {
            centroid[0] += v[0];
            centroid[1] += v[1];
        }
        centroid[0] /= 4.0;
        centroid[1] /= 4.0;

        std::sort(topVertices.begin(), topVertices.end(), [&centroid](const std::vector<double>& a, const std::vector<double>& b) {
            return atan2(a[1] - centroid[1], a[0] - centroid[0]) < atan2(b[1] - centroid[1], b[0] - centroid[0]);
            });

        // 计算面积
        double area = 0.5 * std::abs(
            topVertices[0][0] * (topVertices[1][1] - topVertices[3][1]) +
            topVertices[1][0] * (topVertices[2][1] - topVertices[0][1]) +
            topVertices[2][0] * (topVertices[3][1] - topVertices[1][1]) +
            topVertices[3][0] * (topVertices[0][1] - topVertices[2][1])
        );

        // 计算体积
        double volume = area * (vertices[7][2] - vertices[0][2]); // zmax - zmin
        volumes.push_back(volume);
    }

    return Field(volumes);
}

//计算cell个数
size_t Cells::size() const {
    return cells.size();
}


// Mesh 类构造函数
Mesh::Mesh(const std::string& filename) {
    std::vector<std::vector<double>> nodes0;
    parseNodes(filename, nodes0);
    nodes = nodes0;
    std::vector<std::vector<int>> faceInfo;

    parseFaces(filename, faceInfo);
    Faces Faces0(nodes0, faceInfo); // 在这里赋值
    faces = Faces0;
    std::vector<std::vector<int>> cellsData0;
    std::vector<std::vector<int>> celltreeData0;
    parseFacetree(filename, celltreeData0);
    cellsData0 = processCelldata(faceInfo, celltreeData0);
    Cells Cells0(nodes0, cellsData0); // 在这里赋值
    cells = Cells0;
}

// 计算所有面中心点
Field Mesh::calculateAllfaceCenters() const {
    return faces.calculateAllCenters();
}

// 计算所有面法向量
Field Mesh::calculateAllfaceNormals() const {
    return faces.calculateAllNormals();
}

// 计算所有面面积
Field Mesh::calculateAllfaceAreas() const {
    return faces.calculateAllAreas();
}

// 计算所有单元中心点
Field Mesh::calculateAllcellCenters() const {
    return cells.calculateAllCenters();
}

// 计算所有单元体积
Field Mesh::calculateAllcellVolumes() const {
    return cells.calculateAllVolumes();
}

// 访问索引为 i 的面
const std::vector<int>& Mesh::getFace(size_t i) const {
    return faces.getFace(i);
}
const std::vector<Scalar>& Mesh::getPoint(size_t i) const {
    return nodes[i];
}

// 修改 faces 中zoneid为i的行的边界条件
void Mesh::setBctypeForFace(int zoneIndex, int bctype) {
    faces.setBctypeForFace(zoneIndex, bctype);
}

// 返回面和单元体数量
size_t Mesh::numberOfFaces() const {
    return faces.size();
}

size_t Mesh::numberOfCells() const {
    return cells.size();
}

size_t Mesh::numberOfPoints() const {
    return nodes.size();
}