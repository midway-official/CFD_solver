#include "CFDmath.h"  // 确保包含头文件
#include "FVM.h"

using namespace std;

string filename = "small.msh"; // 输入文件名
// 计算速度场残差和
Scalar computeResidualSumU(const std::vector<std::vector<Scalar>>& u0, const std::vector<std::vector<Scalar>>& u1) {
    if (u0.size() != u1.size()) {
        std::cerr << "Error: u0 and u1 must have the same size." << std::endl;
        return -1.0;
    }
    
    Scalar residualSum = 0.0;
    
    for (size_t i = 0; i < u0.size(); ++i) {
        // 计算每个点上的残差（两个三维向量的欧氏距离）
        Scalar dx = u0[i][0] - u1[i][0];
        Scalar dy = u0[i][1] - u1[i][1];
        Scalar dz = u0[i][2] - u1[i][2];
        residualSum += std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    
    return residualSum;
}

// 计算压力场残差和
Scalar computeResidualSumP(const std::vector<Scalar>& p0, const std::vector<Scalar>& p1) {
    if (p0.size() != p1.size()) {
        std::cerr << "Error: p0 and p1 must have the same size." << std::endl;
        return -1.0;
    }

    Scalar residualSum = 0.0;

    for (size_t i = 0; i < p0.size(); ++i) {
        // 计算每个点上的压力残差
        residualSum += std::abs(p0[i] - p1[i]);
    }

    return residualSum;
}


int main() {
    cout << "开始读取网格" << endl;
    Mesh mesh(filename);
    cout << "完成网格读取, mesh对象创建" << endl;
    cout << "网格基本信息" << endl;
    //计算网格几何信息
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
    Point U_wall(vector<Scalar>{-0.001, 0.0, 0.0});
    
    // 初始猜解初始化
    Field P0(mesh.numberOfCells(), 0.0);
    Field U0(mesh.numberOfCells(), std::vector<Scalar>{0.0, 0.0, 0.0});
    
    // 创建矩阵和源项（仅初始化一次）
    std::vector<std::vector<Scalar>> A0(mesh.numberOfCells(), std::vector<Scalar>(mesh.numberOfCells(), 0.0));
    std::vector<std::vector<Scalar>> A1(mesh.numberOfCells(), std::vector<Scalar>(mesh.numberOfCells(), 0.0));
    std::vector<std::vector<Scalar>> b0(mesh.numberOfCells(), std::vector<Scalar>(3, 0.0));
    std::vector<Scalar> b1(mesh.numberOfCells(),  0.0);
    std::vector<Scalar> scalarVector; // 声明用于存储散度值的向量
for (int iter = 0; iter < 1000; ++iter) {
        // 将 A0 和 A1 重置为全零
        // 重置 A0 和 A1
        for (auto& row : A0) {
            std::fill(row.begin(), row.end(), 0.0); // 将每一行重置为0.0
        }
        for (auto& row : A1) {
            std::fill(row.begin(), row.end(), 0.0); // 将每一行重置为0.0
        }

        // 重置 b0
        for (auto& row : b0) {
            std::fill(row.begin(), row.end(), 0.0); // 将每一行的元素重置为0.0
        }

        // 重置 b1
        std::fill(b1.begin(), b1.end(), 0.0); // 将一维向量 b1 中的所有元素重置为0.0
        // 1. 解动量方程
        Field gradientP = calculateGradient(P0, mesh, 0.0);
        
        Field DvP = (-1 * gradientP * faceVol);
        
        // 组装 b0（只更新源项 b0，而不是重新初始化）
        b0 = DvP.getPointVector();
        
        // 离散对流项, 更新 A0 和 b0
        separateConvective(mesh, 1, A0, b0, U0, U_wall);  // 更新 A0, b0
        separateLaplaceU(mesh,0.0012,A0,b0,U_wall);
        // 使用三维高斯赛德尔求解器解动量方程
        GaussSeidel3D gs3D(0.001,1000,0.8); // 参数收敛精度和最大迭代次数,亚松弛修正
        //获得新场
        auto solution = gs3D.solve(A0, b0);
        Field U1(solution);
        

        // 2. 计算新的速度场的散度
        Field DivU = calculateDivergence(U1, mesh, U_wall);
        b1 = (DivU *(-1)*faceVol).getScalarVector();
        
        // 3. 组装 A1 和 b1（泊松方程，更新压力）
        // 在更新 A1 之前将其重置为全零
        
        separateLaplaceP(mesh, 1, A1, b1, A0);  // 更新 A1, b1
        
        // 4.使用高斯赛德尔求解器解压力修正方程
        GaussSeidel gs(10e-4,1000);
        auto solution2 = gs.solve(A1, b1);  // 注意这里解的是 A1, b1
        //求解压力修正
        Field pc(solution2);
        //求解速度修正
        Field Divpc=calculateGradientPc(pc,mesh);
        Field uc= computeVelocityCorrection(Divpc,A0,faceVol);
        //残差计算
        Scalar resU=computeResidualSumU(U0.getPointVector(),U1.getPointVector());
        Scalar resP=computeResidualSumP(P0.getScalarVector(),(P0+0.01*pc).getScalarVector());
        // 5. 更新压力场和速度场
        P0= P0+0.01*pc;
        U0= U1+uc;
        
        // 计算新的速度场散度以判断是否收敛
        Field divU1 = calculateDivergence(U0, mesh, U_wall);
         // 将散度值存储到向量中
        scalarVector = divU1.getScalarVector();
        Scalar sum = 0.0; // 初始化和为零

        for (const auto& value : scalarVector) {
            sum += value; // 累加每个元素 
        }

        // 打印和
         std::cout << "速度散度场的和(DivU): " << sum <<" 速度场残差(RESU): "<<resU<<" 压力场残差(RESP): "<<resP<< std::endl;
        
        // 判断收敛条件：散度的最大值是否小于阈值
        if (!scalarVector.empty()) {
         

            // 判断收敛条件
            if (abs(sum) < 0.0001 ) {
                std::cout << "SIMPLE算法在第 " << iter + 1 << " 次迭代后收敛。" << std::endl;
               
                break;
            }
        } else {
            std::cerr << "错误：新的速度场散度计算返回了空向量。" << std::endl;
            break;
        }

      
  
    
    }

    U0.print();
    
    return 0;
}
