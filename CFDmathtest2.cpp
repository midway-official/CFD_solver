
#include "CFDmath.h"  // 确保包含头文件
#include "FVM.h"
using namespace std;
string filename = "2.msh"; // 输入文件名

int main(){
     cout << "开始读取网格"  << endl;
    Mesh mesh(filename);
     cout << "完成网格读取,mesh对象创建"  << endl;
     cout << "网格基本信息"  << endl;
    Field facecenter= mesh.calculateAllfaceCenters();
    Field facearea=mesh.calculateAllfaceAreas();
    Field facenorm=mesh.calculateAllfaceNormals();
    Field faceVol=mesh.calculateAllcellVolumes();
    MeshAnalyzer analyzer;

    // 分析 mesh 并打印 zone 信息
 
    //修改出流边界，零压梯度
    mesh.setBctypeForFace(11,36);
     auto zoneInfo = analyzer.analyzeMesh(mesh);
    analyzer.printZoneInfo(zoneInfo);
     cout << "场计算测试"  << endl;
    

    //梯度散度计算测试
    Field P(mesh.numberOfCells(), 0.0);
    Field U(mesh.numberOfCells(), vector<Scalar>{0.0, 0.0, 0.0});
    
    Point U_wall(vector<Scalar>{-1.0, 0.0, 0.0});
    Field gradientP=calculateGradient(P, mesh , 0.0);
    cout << "场梯度计算测试完成"  << endl;
    Field divergenceU=calculateDivergence(U, mesh,U_wall);
    cout << "场散度计算测试完成"  << endl;
    cout << "Number of cells: " << mesh.numberOfCells()  << endl;
    return 0;
    }
        