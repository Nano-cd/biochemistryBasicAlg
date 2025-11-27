好的，遵照您的要求，我为您生成一个详细的 `README.md` 文件。这个文件将全面介绍您提供的C++代码库，解释每个模块的功能、实现逻辑以及如何使用它们。

---

# Biochemical Analysis Calculation Engine in C++

这是一个用 C++ 实现的综合性生化分析计算引擎。该项目旨在模拟现代全自动生化分析仪的核心软件算法，涵盖了从原始吸光度数据处理到浓度计算的完整流程。它为开发、测试和理解生化检测方法背后的数学原理提供了一个强大的基础。

## ✨ 主要功能

- **多种反应类型计算**:
  - **终点法 (Endpoint Method)**: 计算反应完成时的吸光度变化。
  - **固定时间法 (Fixed-Time Method)**: 计算在固定时间间隔内的反应速率。
  - **动力学法 (Kinetic Method)**: 复杂的酶活性计算，包括底物耗尽监测、线性度检查和酶线性扩展。
  - **速率法 (Rate Method)**: 简化的动力学法，使用最小二乘法拟合反应速率。
- **强大的校准模型**:
  - **线性校准 (Linear Calibration)**: 支持单点（K因数）、两点和多点线性回归。
  - **非线性校准 (Non-Linear Calibration)**: 支持多种高级免疫项目拟合模型，包括 `Polynomial 5P`, `Parabola`, `Logistic 4P`, `Exponential 5P` 和 `Spline`（分段线性插值）。
- **完善的质量控制与验证**:
  - **校准验证 (Calibration Validation)**: 全面检查校准曲线的质量，包括重复性(CV%)、空白吸光度、灵敏度、斜率/因子和线性相关性(R²)。
- **模块化设计**:
  - 每个计算模块都封装在独立的类中，职责清晰，易于维护和扩展。
- **纯C++实现**:
  - 仅依赖C++标准库，无第三方依赖，具有良好的跨平台性。包含一个基础的线性代数工具类 (`MatrixUtils`) 以支持非线性拟合。

---

## 🚀 如何编译与运行

该项目是一个单文件 C++ 应用程序，不依赖任何外部库。您可以使用任何支持 C++11 或更高标准的编译器进行编译。

例如，使用 g++:
```bash
g++ -std=c++11 -o analysis_engine your_source_file.cpp
```

然后运行生成的可执行文件:
```bash
./analysis_engine
```

程序将执行 `main` 函数中预设的演示场景，并打印出详细的计算过程和结果。

---

## 📚 模块详解与用法

### 1. 终点法 (`EndpointCalculator`)

- **作用**: 用于测量反应达到终点或平衡时的总吸光度变化。适用于总蛋白、白蛋白等项目。
- **实现逻辑**:
  - 核心公式: `R = Ai - K * Ab`
    - `Ai`: 反应终点（L-M点）的平均吸光度。
    - `Ab`: 空白读数（N-P点）的平均吸光度。
    - `K`: 体积修正因子，用于校正因加入试剂2导致体积变化而引起的稀释效应。
  - 支持样本空白校正: `R' = R_sample - R_sample_blank`。
- **用法**:
  ```cpp
  #include "YourHeader.h" // 包含定义
  
  EndpointCalculator calc;
  Vector raw_curve = { ... }; // 反应曲线数据
  
  // 场景: 双试剂，带体积修正
  EndpointResult result = calc.calculate_raw_response(raw_curve, 5, 15, 50, 60, 0.8);
  std::cout << "最终反应度 R: " << result.R << std::endl;
  
  // 场景: 样本空白校正
  Vector sample_curve = { ... };
  Vector sample_blank_curve = { ... };
  double corrected_R = calc.calculate_corrected_response(sample_curve, sample_blank_curve, 5, 15, 50, 60, 0.8);
  std::cout << "校正后反应度 R': " << corrected_R << std::endl;
  ```

### 2. 固定时间法 (`FixedTimeCalculator`)

- **作用**: 用于测量在特定时间段内反应速率恒定的项目，如肌酐（苦味酸法）。
- **实现逻辑**:
  - 核心公式: `R = 60 * (Slope_Reaction - K * Slope_Blank)`
  - `Slope_Reaction`: 主反应区间 (L-M) 的斜率 `(Abs_M - Abs_L) / Time_Diff`。
  - `Slope_Blank`: 空白区间 (N-P) 的斜率。
  - `K`: 体积修正因子。
  - `60`: 将单位从 `Abs/sec` 转换为 `ΔAbs/min`。
- **用法**:
  ```cpp
  FixedTimeCalculator calc;
  Vector raw_curve = { ... };
  double CYCLE_TIME_SEC = 18.0;
  
  FixedTimeResult result = calc.calculate_response(raw_curve, 5, 16, 17, 30, 0.8, CYCLE_TIME_SEC);
  std::cout << "最终反应度 R: " << result.R << std::endl;
  ```

### 3. 动力学法 (`KineticCalculator`)

- **作用**: 用于酶活性检测，是所有方法中逻辑最复杂的。它能动态适应反应过程中的变化。
- **实现逻辑**:
  1.  **底物耗尽监测**: 从预设区间终点 `M` 向前检查，如果吸光度超出 `depletion_limit`，则自动“读点前移”，缩短计算区间，确保只使用线性部分的数据。
  2.  **速率计算**: 对最终确定的有效区间内的点，使用最小二乘法进行线性回归，计算出最精确的斜率。
  3.  **线性度检查**: 将有效区间分为前后两段，比较两段的斜率差异。如果差异过大，则会报告 `LIN` (Linearity) 报警。
  4.  **酶线性扩展 (EXP)**: 如果主反应区间 `L` 点之后很快就耗尽，导致有效点数不足，算法会自动回溯到延迟期，寻找最快的反应速率作为结果，并报告 `EXP` (Expansion) 标记。
- **用法**:
  ```cpp
  KineticCalculator calc;
  Vector high_activity_curve = { ... };
  double CYCLE_TIME_SEC = 10.0;
  double LIMIT = 0.5; // 吸光度耗尽限值
  int DIR = -1;      // 负反应

  auto result = calc.calculate_response(high_activity_curve, 0, 0, 3, 10, 0, CYCLE_TIME_SEC, LIMIT, DIR);
  std::cout << "反应度: " << result.first << ", 状态: " << result.second << std::endl;
  ```

### 4. 速率法 (`RateMethodCalculator`)

- **作用**: 动力学法的简化版，在指定的固定区间内使用最小二乘法计算反应速率。
- **实现逻辑**: 对 `L` 到 `M` 之间的所有数据点进行线性回归，直接计算斜率。
- **用法**:
  ```cpp
  RateMethodCalculator calc;
  Vector raw_curve = { ... };
  
  RateResult result = calc.calculate_response(raw_curve, 10, 30, 18.0);
  std::cout << "最终反应度 R: " << result.R << " ΔAbs/min" << std::endl;
  
  bool is_linear = calc.linearity_check(raw_curve, 10, 30);
  std::cout << "线性检查: " << (is_linear ? "通过" : "失败") << std::endl;
  ```

### 5. 线性校准 (`LinearCalibrator`)

- **作用**: 将仪器测得的原始反应度 `R` 转换为最终的浓度 `C`。
- **实现逻辑**:
  - **模型**: `C = K * (R - R0)`，其中 `R` 是原始信号值除以一个缩放因子（如10000）得到的吸光度。
  - **单点法**: 用户直接提供 `K` 和 `R0`。
  - **两点法**: 通过两个校准点解出 `K` 和 `R0`。
  - **多点法**: 使用最小二乘法对多个校准点进行线性回归，拟合出最佳的 `K` 和 `R0`。
- **用法**:
  ```cpp
  LinearCalibrator calibrator;

  // 多点校准
  calibrator.fit_multi_point({10.0, 30.0, 50.0}, {1000, 3000, 5000});
  
  // 计算样本浓度
  double concentration = calibrator.calculate_concentration(4000); // 样本原始反应度为4000
  std::cout << "样本浓度: " << concentration << std::endl;
  ```

### 6. 非线性校准 (`NonLinearCalibrator`)

- **作用**: 用于免疫等反应曲线呈 S 型或其他非线性形状的项目。
- **实现逻辑**:
  - `Polynomial 5P`: 将浓度取对数 `ln(C)`，响应值 `R` 做变换后，进行三次多项式拟合。
  - `Parabola`: 将响应值 `R` 拟合为浓度的二次函数 `R = a*C^2 + b*C + R0`。
  - `Logistic 4P`: 经典的四参数逻辑斯蒂模型，通常直接设置由厂家提供的参数。
  - `Spline`: 在相邻的校准点之间进行线性插值，简单有效。
  - **依赖**: 内部使用 `MatrixUtils` 类进行最小二乘法求解。
- **用法**:
  ```cpp
  NonLinearCalibrator calib_nl;

  // 使用 Polynomial 5P 模型拟合
  calib_nl.fit_polynomial_5p({0, 10, 50, 100}, {500, 1200, 2500, 3200});
  
  // 计算样本浓度
  double concentration = calib_nl.calculate_concentration(2000); // 样本响应值为2000
  std::cout << "样本浓度: " << concentration << std::endl;
  ```

### 7. 校准验证 (`CalibrationValidator`)

- **作用**: 在校准完成后，自动检查校准曲线是否满足质量要求。这是确保检测结果可靠性的关键步骤。
- **实现逻辑**:
  - **重复性**: 计算同一点多次测量的 `CV%`，检查是否低于阈值。
  - **空白吸光度**: 检查零浓度点的响应值是否在正常范围内。
  - **灵敏度**: 检查最低浓度点与空白点的响应值差是否足够大。
  - **斜率/因子**: 检查拟合出的线性斜率是否在预设的合理范围内。
  - **线性相关性**: 检查多点校准的 `R²` 值是否足够接近1。
- **用法**:
  ```cpp
  CalibrationValidator validator;
  ValidationParameters params; // 可自定义允收标准
  
  // 模拟一组校准数据 (浓度 -> 响应值列表)
  std::map<double, std::vector<double>> cal_data = {
      {0.0,   {0.051, 0.052}},
      {50.0,  {0.550}},
      {100.0, {1.048}}
  };

  ValidationResult result = validator.validate(cal_data, params);
  std::cout << "验证结果:\n" << result.message << std::endl;
  ```
