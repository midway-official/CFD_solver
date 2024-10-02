
#include "CFDmath.h"  // 确保包含头文件

// 函数声明
void printMaxValue(const vector<vector<int>>& vec) {
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

void printData(const std::vector<std::vector<int>>& faceData) {
    std::cout << "Face Data:" << std::endl;
    size_t maxRows = std::min(faceData.size(), size_t(100000)); // 确保不超过100行
    for (size_t row = 0; row < maxRows; ++row) {
        for (int value : faceData[row]) {
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
void printMaxSubvectorLength(const std::vector<std::vector<int>>& vec) {
    int maxLength = 0;
    for (const auto& subvec : vec) {
        maxLength = std::max(maxLength, static_cast<int>(subvec.size()));
    }
    std::cout<< maxLength<<std::endl;
}
int main() {
    string filename = "2.msh"; // 输入文件名
    vector<vector<double>> nodeData;
    vector<vector<int>> faceData;
    vector<vector<int>> facetreeData;
    parseNodes(filename,nodeData);
    parseFaces(filename,faceData);
    parseFacetree(filename,facetreeData);
    printNodeData(nodeData);
    printData(faceData);
    printMaxValue(faceData);
    // 输出提取的坐标数量
        cout << "总坐标数: " << nodeData.size() << endl;
    // 输出提取的面信息数量
        cout << "总面数: " << faceData.size() << endl;
   
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
    cout << "场加法:\n";
    addedField.print();

    // 测试标量场与点场的乘法
    Field multipliedField = scalarField * pointField; // 标量场与点场的乘法
    cout << "矢量场乘标量场:\n";
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
    //CF.print();
     auto normals = faces.calculateAllNormals();
   Field CN=normals;
    //CN.print();
      auto areas = faces.calculateAllAreas();
   Field a=areas;
    //a.print();
    

    // 使用引用
       std::vector<Point> pointVec2 = { Point({1.0, 2.0, 3.0}), Point({4.0, 5.0, 6.0}) };
    Field field1(pointVec2);
    
    const Point& pa = field1.pointAt(0);  // 引用
    Point pb = field1.pointAt(1);         // 复制
    
    //pb.print();  // 打印 pb

    
    
    auto cellData=processCelldata(faceData,facetreeData);
    printData(cellData);
    printMaxSubvectorLength(cellData);
    std::cout << "Total number of cells: " << cellData.size() << std::endl;
    std::cout << "Total number of points: " << nodeData.size() << std::endl;
     std::cout << "Total number of faces: " << faceData.size() << std::endl;
    
    return 0;
}
