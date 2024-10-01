### CFD数学库功能及用法说明

#### 1. **解析函数**

- **`void parseNodes(const std::string &filename, std::vector<std::vector<double>> &coordinates);`**
  - **用途**: 从给定的.msh文件中读取节点信息并填充`coordinates`向量。
  - **用法**: 
    ```cpp
    std::vector<std::vector<double>> coordinates;
    parseNodes("filename.msh", coordinates);
    ```
  - **说明**: 文件路径为`filename`，`coordinates`将包含读取到的节点坐标。

- **`void parseFaces(const std::string &filename, std::vector<std::vector<int>> &faces);`**
  - **用途**: 从给定的.msh文件中读取面信息并填充`faces`向量。
  - **用法**:
    ```cpp
    std::vector<std::vector<int>> faces;
    parseFaces("filename.msh", faces);
    ```
  - **说明**: 文件路径为`filename`，`faces`将包含读取到的面数据。

#### 2. **类 Point**

- **构造函数**: 
  ```cpp
  Point(const std::vector<double>& coordinates);
  ```
  - **用途**: 初始化一个三维点。
  - **用法**:
    ```cpp
    std::vector<double> coords = {x, y, z};
    Point p(coords);
    ```

- **运算符重载**:
  - **加法**: `Point operator+(const Point& other) const;`
  - **减法**: `Point operator-(const Point& other) const;`
  - **数乘**: `Point operator*(double scalar) const;`
  - **数除**: `Point operator/(double scalar) const;`
  - **叉乘**: `Point operator^(const Point& other) const;`
  - **点乘**: `Scalar dot(const Point& other) const;`

- **其他功能**:
  - **模长**: `double magnitude() const;`
  - **单位化**: `Point normalize() const;`
  - **打印坐标**: `void print() const;`

#### 3. **类 Field**

- **构造函数**:
  - **点场**: `Field(const std::vector<Point>& points);`
  - **标量场**: `Field(const std::vector<Scalar>& scalars);`

- **运算符重载**:
  - **加法**: `Field operator+(const Field& other) const;`
  - **减法**: `Field operator-(const Field& other) const;`
  - **数乘**: `Field operator*(double scalar) const;`
  - **场乘**: `Field operator*(const Field& other) const;`

- **其他功能**:
  - **点乘**: `Field dot(const Field& other) const;`
  - **叉乘**: `Field cross(const Field& other) const;`
  - **模长**: `Field magnitude() const;`
  - **单位化**: `Field normalize() const;`
  - **打印内容**: `void print() const;`

#### 4. **类 Faces**

- **构造函数**:
  ```cpp
  Faces(const std::vector<std::vector<double>>& nodeCoordinates, const std::vector<std::vector<int>>& faceInfo);
  ```
  - **用途**: 初始化面集对象。
  - **用法**:
    ```cpp
    Faces faces(nodeCoordinates, faceInfo);
    ```

- **其他功能**:
  - **打印面信息**: `void print() const;`
  - **计算中心点**: `Field calculateAllCenters() const;`
  - **计算法向量**: `Field calculateAllNormals() const;`
  - **计算面积**: `Field calculateAllAreas() const;`

### 总结
此CFD数学库提供了一系列工具，便于CFD专业人士在数值计算中处理点、场和面的信息。通过解析.msh文件获取数据，利用类和函数进行必要的数学运算，极大地简化了CFD建模与计算的过程。
