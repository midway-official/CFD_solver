#include "CFDmath.h"  // 确保包含头文件
#include "FVM.h"
using namespace std;

string filename = "small.msh"; // 输入文件名

// 打印矩阵 A
void printMatrix(const vector<vector<Scalar>>& A) {
    for (const auto& row : A) {
        for (const auto& value : row) {
            cout << value << "\t"; // 使用制表符分隔
        }
        cout << endl;
    }
}

// 判断矩阵是否为对角主导
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

int main() {
    cout << "开始读取网格" << endl;
    Mesh mesh(filename);
    cout << "完成网格读取, mesh对象创建" << endl;
    cout << "网格基本信息" << endl;

    Field facecenter = mesh.calculateAllfaceCenters();
    Field facearea = mesh.calculateAllfaceAreas();
    Field facenorm = mesh.calculateAllfaceNormals();
    Field faceVol = mesh.calculateAllcellVolumes();
    size_t n = mesh.numberOfCells();
    MeshAnalyzer analyzer;

    // 修改出流边界，零压梯度
    mesh.setBctypeForFace(11, 36);
    auto zoneInfo = analyzer.analyzeMesh(mesh);
    analyzer.printZoneInfo(zoneInfo);
    
    // 边界条件
    Point U_wall(vector<Scalar>{-1.0, 0.0, 0.0});
    
    // 初始猜解初始化
    Field P0(mesh.numberOfCells(), 0.0);
    Field U0(mesh.numberOfCells(), std::vector<Scalar>{0.0, 0.0, 0.0});
    
    // 创建矩阵和源项（仅初始化一次）
    std::vector<std::vector<Scalar>> A0(mesh.numberOfCells(), std::vector<Scalar>(mesh.numberOfCells(), 0.0));
    std::vector<std::vector<Scalar>> A1(mesh.numberOfCells(), std::vector<Scalar>(mesh.numberOfCells(), 0.0));

    std::vector<Scalar> scalarVector; // 声明用于存储散度值的向量
    for (int iter = 0; iter < 100; ++iter) {
        // 将 A0 和 A1 重置为全零
        std::cout << "当前迭代轮数: " << iter + 1 << std::endl; // 打印当前轮数
        std::fill(A0.begin(), A0.end(), std::vector<Scalar>(mesh.numberOfCells(), 0.0));
        std::fill(A1.begin(), A1.end(), std::vector<Scalar>(mesh.numberOfCells(), 0.0));

        // 1. 解动量方程
        Field gradientP = calculateGradient(P0, mesh, 0.0);
        Field DvP = (-1 * gradientP * faceVol);
        
        // 组装 b0（只更新源项 b0，而不是重新初始化）
        std::vector<std::vector<Scalar>> b0 = DvP.getPointVector();
        
        // 离散对流项, 更新 A0 和 b0
        separateConvective(mesh, 1, A0, b0, U0, U_wall);  // 更新 A0, b0
        separateLaplaceU(mesh,0.0012,A0,b0,U_wall);
        // 使用高斯赛德尔求解器解动量方程
        GaussSeidel3D gs3D; // 参数为收敛精度和最大迭代次数
        auto solution = gs3D.solve(A0, b0);
        Field U1(solution);
        U1.print();

        // 2. 计算新的速度场的散度
        Field DivU = calculateDivergence(U1, mesh, U_wall);
        
        // 检查散度是否有效
        if (DivU.getScalarVector().empty()) {
            std::cerr << "错误：散度计算返回了空向量。" << std::endl;
            break;
        }

        // 将散度值存储到向量中
        scalarVector = DivU.getScalarVector();
        Scalar sum = 0.0; // 初始化和为零

        for (const auto& value : scalarVector) {
            sum += value; // 累加每个元素 
        }

        // 打印和
         std::cout << "scalarVector 的和: " << sum << std::endl;
        std::vector<Scalar> b1 = (DivU * -1).getScalarVector();
        
        // 3. 组装 A1 和 b1（泊松方程，更新压力）
        // 在更新 A1 之前将其重置为全零
        std::fill(A1.begin(), A1.end(), std::vector<Scalar>(mesh.numberOfCells(), 0.0));
        separateLaplaceP(mesh, 1, A1, b1, A0);  // 更新 A1, b1
        
        // 使用高斯赛德尔求解器解压力修正方程
        GaussSeidel gs;
        auto solution2 = gs.solve(A1, b1);  // 注意这里解的是 A1, b1
        Field p_(solution2);
        
        // 4. 计算速度修正
        Field divp_ = calculateGradient(p_, mesh, 0.0);
        Field u_ = computeVelocityCorrection(divp_, A0, faceVol);
        
        // 5. 更新压力场和速度场
        Field P1 = P0 + 0.3*p_;
        Field U1_new = U1 + u_;
        
        // 计算新的速度场散度以判断是否收敛
        Field divU1 = calculateDivergence(U1_new, mesh, U_wall);
        
        // 判断收敛条件：散度的最大值是否小于阈值
        if (!scalarVector.empty()) {
         

            // 判断收敛条件
            if (abs(sum) < 0.1 ) {
                std::cout << "SIMPLE算法在第 " << iter + 1 << " 次迭代后收敛。" << std::endl;
               
                break;
            }
        } else {
            std::cerr << "错误：新的速度场散度计算返回了空向量。" << std::endl;
            break;
        }

        // 更新初始猜解以进行下一轮迭代
        P0 = P1;
        U0 = U1_new;
  
    
    }

    U0.print();
    
    return 0;
}
