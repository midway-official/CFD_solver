
#include "CFDmath.h"  // 确保包含头文件
#include "FVM.h"
using namespace std;
string filename = "2.msh"; // 输入文件名

Scalar findMinValue(const Field& field) {
   

    // 假设标量场非空，初始化最小值为第一个元素
    Scalar minValue = field.scalarAt(0);
     
    // 遍历标量场，查找最小值
    
    for (size_t i = 1; i < field.size(); ++i) {
        if (field.scalarAt(i) < minValue) {
            minValue = field.scalarAt(i);
        }
    }

    return minValue;
}


int main(){
    Mesh mesh(filename);
     cout << "完成网格读取,mesh对象创建"  << endl;
    Field facecenter= mesh.calculateAllfaceCenters();
    Field facearea=mesh.calculateAllfaceAreas();
    Field facenorm=mesh.calculateAllfaceNormals();
    Field faceVol=mesh.calculateAllcellVolumes();
     cout << "场计算测试"  << endl;
     double minv =findMinValue(faceVol);
     Field cp =facenorm/facearea;
    Field cd =facearea/facearea;
    cout << "cell最小体积: " << minv  << endl;

    Field P(mesh.numberOfCells(), 0.0);
    Field U(mesh.numberOfCells(), vector<Scalar>{0.0, 0.0, 0.0});

    
    Field gradientP=calculateGradient(P, mesh);
    cout << "场梯度计算测试完成"  << endl;
     //gradientP.print();
    Field divergenceU=calculateDivergence(U, mesh);
    cout << "场散度计算测试完成"  << endl;
    cout << "Number of cells: " << mesh.numberOfCells()  << endl;
    return 0;
    }
        