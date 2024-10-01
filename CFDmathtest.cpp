
#include "CFDmath.h"  // 确保包含头文件

// 函数声明
void printMaxValue(const std::vector<std::vector<int>>& vec) {
    if (vec.empty() || vec[0].empty()) {
        std::cout << "The vector is empty." << std::endl;
        return;
    }

    int maxValue = std::numeric_limits<int>::min(); // 初始化为最小整数

    for (const auto& innerVec : vec) {
        for (int value : innerVec) {
            if (value > maxValue) {
                maxValue = value;
            }
        }
    }

    std::cout << "Maximum value: " << maxValue << std::endl;
}
void printNodeData(const std::vector<std::vector<double>>& nodeData) {
    std::cout << "Node Data:" << std::endl;
    for (const auto& node : nodeData) {
        for (double value : node) {
            std::cout << value << " ";
        }
        std::cout << std::endl; // 换行
    }
}

void printFaceData(const std::vector<std::vector<int>>& faceData) {
    std::cout << "Face Data:" << std::endl;
    for (const auto& face : faceData) {
        for (int value : face) {
            std::cout << value << " ";
        }
        std::cout << std::endl; // 换行
    }
}
void printFaces(const std::vector<std::vector<int>>& faces) {
    for (const auto& face : faces) {
        std::cout << "(";
        for (size_t i = 0; i < face.size(); ++i) {
            std::cout << face[i];
            if (i < face.size() - 1) {
                std::cout << ", "; // 添加逗号分隔
            }
        }
        std::cout << ")" << std::endl; // 每个子向量换行
    }
}
int main() {
    std::string filename = "2.msh"; // 输入文件名
    std::vector<std::vector<double>> nodeData;
    std::vector<std::vector<int>> faceData;
    
    parseNodes(filename,nodeData);
    parseFaces(filename,faceData);
    printNodeData(nodeData);
    printMaxValue(faceData);
    // 输出提取的坐标数量
        std::cout << "Total number of coordinates: " << nodeData.size() << std::endl;
    // 输出提取的面信息数量
        std::cout << "Total number of faces: " << faceData.size() << std::endl;
   
    Point p1({1.0, 2.0, 3.0});
    Point p2({4.0, 5.0, 6.0});

    Point p3 = p1 + p2;
    Point p4 = p1 - p2;
    Point p5 = p1 * 2.0;
    Point p6 = p1 / 2.0;
    Point p7 = p1 ^ p2;
    double mag = p1.magnitude();
    Point unit = p1.normalize();

    std::cout << "p1: "; p1.print();
    std::cout << "p2: "; p2.print();
    std::cout << "p3 (p1 + p2): "; p3.print();
    std::cout << "p4 (p1 - p2): "; p4.print();
    std::cout << "p5 (p1 * 2): "; p5.print();
    std::cout << "p6 (p1 / 2): "; p6.print();
    std::cout << "p7 (p1 ^ p2): "; p7.print();
    std::cout << "Magnitude of p1: " << mag << "\n";
    std::cout << "Unit vector of p1: "; unit.print();

    Faces faces(nodeData, faceData);
        //faces.printFaces();
    // 计算并打印所有面的中心点
    auto centers = faces.calculateAllCenters();
   /* for (const auto& center : centers) {
        std::cout << "Center: ";
        center.print();
    }*/
    
    return 0;
}
