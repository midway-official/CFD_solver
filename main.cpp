
#include "CFDmath.h"  // 确保包含头文件
#include "FVM.h"
using namespace std;
string filename = "2.msh"; // 输入文件名
// 打印矩阵 A
void printMatrix(const vector<vector<Scalar>>& A) {
    for (const auto& row : A) {
        for (const auto& value : row) {
            cout << value << "\t"; // 使用制表符分隔
        }
        cout << endl;
    }
}
bool isDiagonallyDominant(const std::vector<std::vector<double>>& A) {
    int n = A.size(); // 矩阵的行数
    for (int i = 0; i < n; ++i) {
        double diagonalElement = fabs(A[i][i]); // 对角线元素的绝对值
        double sumOfRow = 0.0; // 当前行的其他元素绝对值之和

        // 计算当前行的其他元素的绝对值之和
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                sumOfRow += fabs(A[i][j]);
            }
        }

        // 检查对角主导性
        if (diagonalElement <= sumOfRow) {
            return false; // 如果某一行不满足对角主导性，返回 false
        }
    }
    return true; // 如果所有行都满足对角主导性，返回 true
}
// 打印向量 b
void printVector(const vector<Scalar>& b) {
    for (const auto& value : b) {
        cout << value << "\t"; // 使用制表符分隔
    }
    cout << endl;
}
int countZeroDiagonal(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();  // 获取矩阵的大小
    int count = 0;

    // 遍历对角线元素
    for (int i = 0; i < n; ++i) {
        if (matrix[i][i] == 0.0) {
            ++count;  // 统计对角元素为零的个数
        }
    }

    return count;
}
int main(){
    /*
     cout << "开始读取网格"  << endl;
    Mesh mesh(filename);
     cout << "完成网格读取,mesh对象创建"  << endl;
     cout << "网格基本信息"  << endl;
    Field facecenter= mesh.calculateAllfaceCenters();
    Field facearea=mesh.calculateAllfaceAreas();
    Field facenorm=mesh.calculateAllfaceNormals();
    Field faceVol=mesh.calculateAllcellVolumes();
    size_t n = mesh.numberOfCells();
    MeshAnalyzer analyzer;

    // 分析 mesh 并打印 zone 信息
 
    //修改出流边界，零压梯度
    mesh.setBctypeForFace(11,36);
     auto zoneInfo = analyzer.analyzeMesh(mesh);
    analyzer.printZoneInfo(zoneInfo);
    //边界条件
    Point U_wall(vector<Scalar>{-1.0, 0.0, 0.0});
    

   
   //simple算法
   //初始化
    Field P0(mesh.numberOfCells(), 0.0);
    Field U0(mesh.numberOfCells(), vector<Scalar>{-1, 0.0, 0.0});    
    vector<vector<Scalar>> A0(n, std::vector<Scalar>(n, 0.0));
    vector<vector<Scalar>> A1(n, std::vector<Scalar>(n, 0.0));
    vector<Scalar> b0(n, 0);
    vector<Scalar> b1(n, 0);

    //解动量方程
    //求压力梯度
    Field gradientP=calculateGradient(P0, mesh , 0.0);
    Field DvP=(-1*gradientP*faceVol); // 确保 b 的大小与单元体数量相同
    //组装b
 
    //离散对流项,更新A0和b0
    separateConvective(mesh,1,A0,b0,U0,U_wall);
    vector<Scalar> b0X=addVectors(DvP.getXComponent(),b0);
    vector<Scalar> b0Y=addVectors(DvP.getYComponent(),b0);
    vector<Scalar> b0Z=addVectors(DvP.getZComponent(),b0);
    GaussSeidel solver(10e-5,10);
    
    
    GaussSeidel3D gs(0.1,20,0.1); // 创建高斯-赛德尔3D求解器实例
    auto solution = gs.solve(A0,b0X,b0Y,b0Z); // 调用求解器

    // 输出结果
    std::cout << "\nSolution:" << std::endl;
    for (size_t i = 0; i < solution.size(); i++) {
        std::cout << "x[" << i << "] = (" 
                  << solution[i][0] << ", " << solution[i][1] << ", " << solution[i][2] << ")" << std::endl; // 输出解向量
    }
    */
     // 示例三维矢量方程组
    int n = 2; // 设定有两个方程
    std::vector<std::vector<double>> A = {
        {2, 1}, // 对应于 A[0]
        {1, 2}  // 对应于 A[1]
    };

    // 使用三个分量表示 b 的三维矢量
    std::vector<double> b0X = {5, 4}; // 常数矢量分量 X
    std::vector<double> b0Y = {3, 2}; // 常数矢量分量 Y
    std::vector<double> b0Z = {0, 0}; // 常数矢量分量 Z

    GaussSeidel3D gs(1e-10, 1000, 1.0); // 创建高斯-赛德尔3D求解器实例，松弛因子为1.0（无松弛）
    auto solution = gs.solve(A, b0X, b0Y, b0Z); // 调用求解器

    // 输出解中模长最大的向量
    std::vector<double> maxVector(3, 0.0);
    double maxLength = 0.0;

    for (const auto& vec : solution) {
        double length = sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2)); // 计算模长
        if (length > maxLength) {
            maxLength = length;
            maxVector = vec; // 更新最大模长的向量
        }
    }

    // 输出结果
    std::cout << "\nSolution:" << std::endl;
    for (size_t i = 0; i < solution.size(); i++) {
        std::cout << "x[" << i << "] = (" 
                  << solution[i][0] << ", " << solution[i][1] << ", " << solution[i][2] << ")" << std::endl; // 输出解向量
    }

    // 输出模长最大的向量
    std::cout << "\nThe vector with the maximum length is: (" 
              << maxVector[0] << ", " << maxVector[1] << ", " << maxVector[2] << ")" 
              << " with length: " << maxLength << std::endl;

 

    return 0;
    }
        