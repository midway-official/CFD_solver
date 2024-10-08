#include "CFDmath.h"  // 数学库头文件
#include "FVM.h"
using namespace std;

string getBoundaryType(int bctype) {
    map<int, string> boundaryTypes = {
        {0x2, "interior"},
        {0x3, "wall"},
        {0x4, "pressure-inlet"},
        {0x5, "pressure-outlet"},
        {0x7, "symmetry"},
        {0x8, "periodic-shadow"},
        {0x9, "pressure-far-field"},
        {0xa, "velocity-inlet"},
        {0xc, "periodic"},
        {0xe, "fan, porous-jump, radiator"},
        {0x14, "mass-flow-inlet"},
        {0x18, "Interface"},
        {0x1f, "parent"},
        {0x24, "outflow"},
        {0x25, "axis"}
    };

    if (boundaryTypes.find(bctype) != boundaryTypes.end()) {
        return boundaryTypes[bctype];
    } else {
        return "unknown";
    }
}

vector<BoundaryCondition> parseBoundaryConditions(const string& filename) {
    ifstream file(filename);
    string line;
    vector<BoundaryCondition> boundaryConditions;

    // 检查文件是否正常打开
    if (!file.is_open()) {
        cerr << "无法打开文件: " << filename << endl;
        return boundaryConditions;
    }

    // 正则表达式匹配 (39 (zoneid bctype regionName)
    regex pattern(R"\(\d+\s+\((\d+)\s+(\w+)\s+(\w+)\)\s+\(\))");
    smatch matches;

    // 逐行解析文件
    while (getline(file, line)) {
        // 使用正则表达式查找符合格式的行
        if (regex_search(line, matches, pattern)) {
            BoundaryCondition bc;
            bc.zoneid = stoi(matches[1].str()); // zoneid
            bc.bctype = stoi(matches[2].str()); // bctype (十六进制形式)
            bc.regionName = matches[3].str();    // regionName

            boundaryConditions.push_back(bc);
        }
    }

    return boundaryConditions;
}

Field calculateGradient(Field& P, Mesh& mesh) {
    Field gradientP(mesh.numberOfCells(), vector<double>{0, 0, 0});
    Field facearea = mesh.calculateAllfaceAreas();
    Field facenorm = mesh.calculateAllfaceNormals();
    Field cellVol = mesh.calculateAllcellVolumes();


    for (size_t i = 0; i < mesh.numberOfFaces(); ++i) {
       

        if (mesh.getFace(i)[1] != 2) {
          
            continue; // 如果不是2（内部面），跳过当前迭代
        }

        size_t c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
        size_t c1 = mesh.getFace(i)[7] - 1;
        
       
        // 确保 c0 和 c1 在有效范围内
        if (c0 >= P.size() || c1 >= P.size()) {
            std::cerr << "Error: Cell index out of range. c0: " << c0 << ", c1: " << c1 << ", P size: " << P.size() << std::endl;
            throw out_of_range("Cell超出索引.");
        }
        
        Scalar Pf = -(P.scalarAt(c0) + P.scalarAt(c1)) / 2.0; // 法向量指向c0，故加负号表示朝外
        double A = facearea.scalarAt(i); // 当前面面积
        

        Point gradientPc0 = gradientP.pointAt(c0);
        Point gradientPc1 = gradientP.pointAt(c1);
        
        
        gradientP.setPoint(c0, gradientPc0 + facenorm.pointAt(i) * Pf * A);
        gradientP.setPoint(c1, gradientPc1 + facenorm.pointAt(i) * Pf * A * (-1.0)); // c0 c1外向面法向量相反
    }

    std::cout << "Final gradient calculated." << std::endl;

    return gradientP / cellVol;
}


Field calculateDivergence(Field& U, Mesh& mesh) {
    Field divergenceU(mesh.numberOfCells(),0.0);
    Field facearea = mesh.calculateAllfaceAreas();
    Field facenorm = mesh.calculateAllfaceNormals();
    Field cellVol = mesh.calculateAllcellVolumes();


     for (size_t i = 0; i < mesh.numberOfFaces(); ++i) {
       
            
        if (mesh.getFace(i)[1] != 2) {
          
            continue; // 如果不是2（内部面），跳过当前迭代
        }

        size_t c0 = mesh.getFace(i)[6] - 1; // 编号减1转换为索引
        size_t c1 = mesh.getFace(i)[7] - 1;
        
       
       
        
       Point Uf = (U.pointAt(c0) + U.pointAt(c1))*(-0.5) ; // 计算面上的速度，加面法向量指入c0负号


        Scalar A = facearea.scalarAt(i); // 当前面面积
      
        // 散度更新
        Scalar divergenceContributionC0 = A * (Uf.dot(facenorm.pointAt(i)));
        Scalar divergenceContributionC1 = A * (Uf.dot(facenorm.pointAt(i)));

        

        divergenceU.setScalar(c0, divergenceU.scalarAt(c0) + divergenceContributionC0);
        divergenceU.setScalar(c1, divergenceU.scalarAt(c1) - divergenceContributionC1); // 法向量相反
    }

   
    std::cout << "Final divergence calculated." << std::endl;

    

    return divergenceU / cellVol;
}