#include "CFDmath.h"
#include <iostream>

int main() {
    try {
        // 测试 n_head_to_array
        std::string n_head_input = "1 2 3 4 5 0x1A";
        std::vector<int> n_result = n_head_to_array(n_head_input);
        std::cout << "n_head_to_array 结果: ";
        for (int val : n_result) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // 测试 f_head_to_array
        std::string f_head_input = "1 23 23c2 234 0x5";
        std::vector<int> f_result = f_head_to_array(f_head_input);
        std::cout << "f_head_to_array 结果: ";
        for (int val : f_result) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        // 测试 nodes_coordinate_to_array
        std::string nodes_input = "1.0 2.0 3.0 4.0 5.0 6.0";
        auto nodes_result = nodes_coordinate_to_array(nodes_input, 3);
        std::cout << "nodes_coordinate_to_array 结果: " << std::endl;
        for (const auto& coords : nodes_result) {
            for (const auto& coord : coords) {
                std::cout << coord << " ";
            }
            std::cout << std::endl;
        }

        // 测试 Mesh 类
        Mesh mesh;
       
        mesh.parse_msh_file("2.msh");  // 读取文件
        print2DArraySize(mesh.get_points());
        print2DArraySize(mesh.get_faces());

        //测试nodes类
        Nodes nodes(mesh.get_points()) ;
        nodes.print();
        //测试faces类
        Faces faces(mesh.get_faces(),nodes);
        faces.calculateAllNormals(nodes);
        faces.calculateAllCenters(nodes);
        faces.calculateAllAreas(nodes);
   
    } catch (const std::exception& e) {
        std::cerr << "发生异常: " << e.what() << std::endl;
    }

    return 0;
}
