
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
   
    // 创建一些点
    Point p1({1.0, 2.0, 3.0});
    Point p2({4.0, 5.0, 6.0});
    Point p3({7.0, 8.0, 9.0});
    
    // 创建点场
    std::vector<Point> pointVec = {p1, p2, p3};
    Field pointField(pointVec);

    // 创建标量场
    std::vector<Scalar> scalarVec = {2.0, 3.0, 4.0};
    Field scalarField(scalarVec);

    // 测试加法
    Field addedField = pointField + pointField; // 点场加法
    std::cout << "Added Point Field:\n";
    addedField.print();

    // 测试标量场与点场的乘法
    Field multipliedField = scalarField * pointField; // 标量场与点场的乘法
    std::cout << "Multiplied Point Field:\n";
    multipliedField.print();

    // 测试点场的模长
    Field magnitudesField = pointField.magnitude();
    std::cout << "Magnitudes of Point Field:\n";
    magnitudesField.print();

    // 测试点场的标准化
    Field normalizedField = pointField.normalize();
    std::cout << "Normalized Point Field:\n";
    normalizedField.print();

    // 测试点场之间的点乘
    Field dotField = pointField.dot(pointField); // 自点乘
    std::cout << "Dot Product of Point Fields:\n";
    dotField.print();

    // 测试点场之间的叉乘
    Field crossField = pointField.cross(pointField); // 自叉乘
    std::cout << "Cross Product of Point Fields:\n";
    crossField.print();

    
    Faces faces(nodeData, faceData);
        //faces.printFaces();
    // 计算并打印所有面的中心点
    auto centers = faces.calculateAllCenters();
   Field CF=centers;
    CF.print();
     auto normals = faces.calculateAllNormals();
   Field CN=normals;
    CN.print();
      auto areas = faces.calculateAllAreas();
   Field a=areas;
    a.print();
    
    return 0;
}
