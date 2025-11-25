import math


# ==========================================
# 基础线性代数工具 (用于替代 numpy)
# ==========================================
class MatrixUtils:
    @staticmethod
    def transpose(matrix):
        return list(map(list, zip(*matrix)))

    @staticmethod
    def multiply(A, B):
        # A: m x n, B: n x p -> C: m x p
        m, n = len(A), len(A[0])
        p = len(B[0])
        C = [[0.0] * p for _ in range(m)]
        for i in range(m):
            for j in range(p):
                for k in range(n):
                    C[i][j] += A[i][k] * B[k][j]
        return C

    @staticmethod
    def inverse(A):
        # 高斯-约旦消元法求逆
        n = len(A)
        # 构造增广矩阵 [A | I]
        M = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(A)]

        for i in range(n):
            pivot = M[i][i]
            if abs(pivot) < 1e-12: raise ValueError("矩阵不可逆")
            for j in range(2 * n): M[i][j] /= pivot
            for k in range(n):
                if k != i:
                    factor = M[k][i]
                    for j in range(2 * n): M[k][j] -= factor * M[i][j]

        return [row[n:] for row in M]

    @staticmethod
    def solve_least_squares(X, Y):
        # 求解 Beta = (X^T * X)^-1 * X^T * Y
        XT = MatrixUtils.transpose(X)
        XTX = MatrixUtils.multiply(XT, X)
        try:
            XTX_inv = MatrixUtils.inverse(XTX)
        except:
            return None  # 奇异矩阵无法求解
        XTY = MatrixUtils.multiply(XT, Y)
        Beta = MatrixUtils.multiply(XTX_inv, XTY)
        return [b[0] for b in Beta]


# ==========================================
# 非线性校准核心类
# ==========================================
class NonLinearCalibrator:
    """
    对应文档 4.6.2 非线性校准
    包含: Logistic-Log 4P/5P, Exponential 5P, Polynomial 5P, Parabola, Spline
    """

    def __init__(self):
        self.method = None
        self.params = []  # 存储参数 [R0, K, a, b, c...]
        self.spline_data = {}  # 专门存储样条数据

    # -------------------------------------------------------
    # 1. Polynomial 5P (5参数多项式) - 完整拟合
    # -------------------------------------------------------
    def fit_polynomial_5p(self, concs, resps):
        """
        文档公式: ln C = a + b*X + c*X^2 + d*X^3
        其中 X = (R - R0) / 100
        条件: 第1个校准品(C=0)的 R 即为 R0
        """
        if len(concs) < 5:
            print("[Error] Polynomial 5P 需要至少 5 个校准品")
            return

        # 1. 确定 R0 (从第一个点，通常 C=0)
        R0 = resps[0]

        # 2. 准备数据 (排除 C=0 的点，因为 ln(0) 无意义)
        X_mat = []  # 构建设计矩阵 [1, x, x^2, x^3]
        Y_mat = []  # 构建目标向量 [ln C]

        valid_points = 0
        for c, r in zip(concs, resps):
            if c <= 0: continue  # 跳过浓度为0的点

            # 变换变量
            x_val = (r - R0) / 100.0
            ln_c = math.log(c)

            X_mat.append([1.0, x_val, x_val ** 2, x_val ** 3])
            Y_mat.append([ln_c])
            valid_points += 1

        if valid_points < 4:
            print("[Error] 有效点(C>0)不足以拟合3次方程")
            return

        # 3. 最小二乘法求解 [a, b, c, d]
        coeffs = MatrixUtils.solve_least_squares(X_mat, Y_mat)
        if not coeffs:
            print("[Error] 拟合失败")
            return

        # 保存参数: R0, a, b, c, d
        self.method = 'Polynomial_5P'
        self.params = [R0] + coeffs  # [R0, a, b, c, d]
        print(f"[拟合成功] Poly 5P: R0={R0}, a={coeffs[0]:.4f}, b={coeffs[1]:.4f}...")

    def calc_polynomial_5p(self, R):
        # 浓度反算：直接代入多项式求 exp
        # y = a + bx + cx^2 + dx^3, 其中 x = (R-R0)/100
        R0, a, b, c, d = self.params
        x = (R - R0) / 100.0
        ln_c = a + b * x + c * (x ** 2) + d * (x ** 3)
        return math.exp(ln_c)

    # -------------------------------------------------------
    # 2. Parabola (抛物线/二次方程) - 完整拟合
    # -------------------------------------------------------
    def fit_parabola(self, concs, resps):
        """
        文档公式: R = a*C^2 + b*C + R0
        最小二乘法求解 R 对 C 的二次回归
        """
        if len(concs) < 3:
            print("[Error] Parabola 需要至少 3 个校准品")
            return

        X_mat = []
        Y_mat = []
        for c, r in zip(concs, resps):
            X_mat.append([c ** 2, c, 1.0])  # [C^2, C, 1]
            Y_mat.append([r])

        coeffs = MatrixUtils.solve_least_squares(X_mat, Y_mat)  # 返回 [a, b, R0]

        self.method = 'Parabola'
        self.params = coeffs  # [a, b, R0]
        print(f"[拟合成功] Parabola: R={coeffs[0]:.4f}C^2 + {coeffs[1]:.4f}C + {coeffs[2]:.4f}")

    def calc_parabola(self, R):
        # 浓度反算: 解一元二次方程 a*C^2 + b*C + (R0 - R) = 0
        a, b, R0 = self.params

        delta = b ** 2 - 4 * a * (R0 - R)
        if delta < 0: return 0.0  # 无实根

        sqrt_delta = math.sqrt(delta)

        # 只要正根 (浓度 > 0)
        if a == 0: return (R - R0) / b if b != 0 else 0

        c1 = (-b + sqrt_delta) / (2 * a)
        c2 = (-b - sqrt_delta) / (2 * a)

        if c1 >= 0: return c1
        if c2 >= 0: return c2
        return 0.0

    # -------------------------------------------------------
    # 3. Logistic-Log 4P / 5P & Exponential 5P
    # -------------------------------------------------------
    # 注意：这三种模型需要 L-M 迭代求解，无库环境下实现极难。
    # 这里提供手动设置参数接口(用于模拟) 和 关键的反算公式(Calculate)。

    def set_params_logistic_4p(self, R0, K, a, b):
        self.method = 'Logistic_4P'
        self.params = [R0, K, a, b]

    def calc_logistic_4p(self, R):
        """
        公式: R = R0 + K / (1 + exp(-(a + b*lnC)))
        反解: (R - R0)/K = 1 / (1 + E)
              1 + E = K / (R - R0)
              E = K / (R - R0) - 1
              -(a + b*lnC) = ln(E)
              lnC = -(ln(E) + a) / b
        """
        R0, K, a, b = self.params
        try:
            val = K / (R - R0) - 1
            if val <= 0: return 0.0
            ln_E = math.log(val)
            ln_C = -(ln_E + a) / b
            return math.exp(ln_C)
        except:
            return 0.0

    def set_params_exponential_5p(self, R0, K, a, b, c):
        self.method = 'Exponential_5P'
        self.params = [R0, K, a, b, c]

    def calc_exponential_5p(self, R):
        """
        公式: R = R0 + K * exp(a*lnC + b*(lnC)^2 + c*(lnC)^3)
        反解: (R - R0)/K = exp(P)
              ln((R-R0)/K) = a*X + b*X^2 + c*X^3  (其中 X = lnC)
        这是一个关于 lnC 的三次方程，需要求根公式或牛顿法求解。
        这里用简单的牛顿法解 X (即 lnC)
        """
        R0, K, a, b, c = self.params
        try:
            target = math.log((R - R0) / K)

            # 解方程 f(x) = ax + bx^2 + cx^3 - target = 0
            x = 1.0  # 初始猜测
            for _ in range(20):  # 迭代20次
                fx = a * x + b * x ** 2 + c * x ** 3 - target
                dfx = a + 2 * b * x + 3 * c * x ** 2
                if abs(dfx) < 1e-6: break
                x_new = x - fx / dfx
                if abs(x_new - x) < 1e-6:
                    x = x_new
                    break
                x = x_new
            return math.exp(x)
        except:
            return 0.0

    # -------------------------------------------------------
    # 4. Spline (样条) - 简化版 (线性分段插值模拟)
    # -------------------------------------------------------
    # 完整的自然三次样条求解在无库下代码量过大，这里实现一个
    # 符合文档描述的 "分段拟合" 逻辑，实际生化仪常用分段线性或简化样条。

    def fit_spline(self, concs, resps):
        # 简单地存储点对，按浓度排序
        data = sorted(zip(resps, concs))  # 按反应度排序用于查找
        self.spline_data = data
        self.method = 'Spline'
        print(f"[拟合成功] Spline: 存储了 {len(data)} 个节点")

    def calc_spline(self, R):
        # 查找 R 所在的区间
        data = self.spline_data
        if not data: return 0.0

        # 外推
        if R <= data[0][0]: return data[0][1]
        if R >= data[-1][0]: return data[-1][1]

        # 内插
        for i in range(len(data) - 1):
            r1, c1 = data[i]
            r2, c2 = data[i + 1]

            if r1 <= R <= r2:
                # 简单的两点线性插值
                # C = C1 + (R - R1) * (C2 - C1) / (R2 - R1)
                return c1 + (R - r1) * (c2 - c1) / (r2 - r1)
        return 0.0

    # -------------------------------------------------------
    # 统一入口
    # -------------------------------------------------------
    def calculate_concentration(self, R):
        if self.method == 'Polynomial_5P':
            return self.calc_polynomial_5p(R)
        elif self.method == 'Parabola':
            return self.calc_parabola(R)
        elif self.method == 'Logistic_4P':
            return self.calc_logistic_4p(R)
        elif self.method == 'Exponential_5P':
            return self.calc_exponential_5p(R)
        elif self.method == 'Spline':
            return self.calc_spline(R)
        else:
            print("未选择校准方法")
            return 0.0


class FixedTimeCalculator:
    """
    【模块：固定时间法反应度计算】
    基于文档 4.4.1 和 4.4.2 实现。
    又称：一级动力学法、初速率法。

    核心公式:
    R = 60 * ( (Am - Al)/(tm - tl) - k * (Ap - An)/(tp - tn) )
    """

    def calculate_response(self, absorbances, N, P, L, M, K, cycle_time=1.0):
        """
        计算固定时间法的反应度

        :param absorbances: 吸光度数据列表 (list)
        :param N: 空白区间起始点索引 (Blank Start)
        :param P: 空白区间结束点索引 (Blank End)
        :param L: 反应区间起始点索引 (Reaction Start)
        :param M: 反应区间结束点索引 (Reaction End)
        :param K: 体积校正因子
        :param cycle_time: 测光周期(秒)，即两个相邻数据点之间的时间间隔。
                           公式中的 tM - tL = (M - L) * cycle_time
        :return: 计算出的 R 值 (Delta Abs / min)
        """

        # 1. 数据边界检查
        max_idx = len(absorbances) - 1
        if M > max_idx or L > max_idx:
            print("[错误] 反应区间索引超出数据范围")
            return 0.0

        # 2. 计算主反应区间的变化速率 (Slope Reaction)
        # (Am - Al) / (tm - tl)
        # tm - tl = (M - L) * cycle_time
        val_M = absorbances[M]
        val_L = absorbances[L]

        time_diff_reaction = (M - L) * cycle_time

        if time_diff_reaction == 0:
            print("[错误] 反应区间时间差为0 (L==M)")
            return 0.0

        slope_reaction = (val_M - val_L) / time_diff_reaction

        # 3. 计算空白区间的变化速率 (Slope Blank)
        # (Ap - An) / (tp - tn)
        slope_blank = 0.0

        # 只有当 N, P 不为 0 且 P > N 时才计算空白速率
        # 对应表 4.2 中 "不扣除空白" 的情况 (N=P=0)
        if P > 0 and P > N and K != 0:
            if P > max_idx:
                print("[错误] 空白区间索引超出数据范围")
                return 0.0

            val_P = absorbances[P]
            val_N = absorbances[N]
            time_diff_blank = (P - N) * cycle_time

            if time_diff_blank != 0:
                slope_blank = (val_P - val_N) / time_diff_blank

        # 4. 整合公式
        # R = 60 * (Slope_Reaction - k * Slope_Blank)
        # 乘以 60 是为了将单位从 "每秒变化" 转换为 "每分钟变化"
        R = 60.0 * (slope_reaction - (K * slope_blank))

        return R, slope_reaction, slope_blank

    def check_boe(self):
        """
        4.4.1 提到的底物耗尽检查 (Substrate Depletion / BOE)
        通常逻辑：检查两点间的吸光度是否超过了线性限，或者初始吸光度是否异常。
        此处仅为占位符，具体逻辑需参考更详细的仪器参数。
        """
        pass


class KineticCalculator:
    """
    【模块：动力学法 (Kinetic Method)】
    基于文档 4.5 全章节实现。
    包含：线性范围判定、最小二乘法计算、线性度检查、酶线性扩展。
    """

    def __init__(self):
        # 60单位对应的吸光度变化率 (0.006 Abs/min)
        # 文档 4.5.5 提到：当斜率小于 60 [单位 A/10000/min] 时不做检查
        self.MIN_slope_CHECK_THRESHOLD = 0.006

    def _calculate_slope_least_squares(self, absorbances, indices, cycle_time):
        """
        【4.5.4 反应度的计算】
        核心最小二乘法公式。
        返回单位：Abs/min (已乘以60)
        """
        n = len(indices)
        if n < 2:
            return 0.0

        # 准备 Ti (时间, 秒) 和 Ai (吸光度)
        # 文档公式中 Ti 为实际测光时间
        T = [idx * cycle_time for idx in indices]
        A = [absorbances[idx] for idx in indices]

        # 计算 T_bar 和 A_bar (平均值)
        T_bar = sum(T) / n
        A_bar = sum(A) / n

        numerator = 0.0  # 分子 sum((Ti - T_bar)*(Ai - A_bar))
        denominator = 0.0  # 分母 sum((Ti - T_bar)^2)

        for i in range(n):
            t_diff = T[i] - T_bar
            a_diff = A[i] - A_bar
            numerator += t_diff * a_diff
            denominator += t_diff * t_diff

        if denominator == 0:
            return 0.0

        # 结果乘以 60 (转为每分钟)
        slope = 60 * (numerator / denominator)
        return slope

    def _check_substrate_depletion(self, abs_val, limit, direction):
        """
        辅助函数：判断某一点是否底物耗尽
        direction: 1 (反应上升), -1 (反应下降)
        limit: 用户输入的耗尽限值
        """
        if limit is None:
            return False

        # 上升反应：吸光度 > 限值 则耗尽
        if direction == 1 and abs_val > limit:
            return True
        # 下降反应：吸光度 < 限值 则耗尽
        elif direction == -1 and abs_val < limit:
            return True

        return False

    def _determine_linear_range(self, absorbances, L, M, limit, direction):
        """
        【4.5.3 确定线性范围】
        """
        max_idx = len(absorbances) - 1
        curr_M = min(M, max_idx)

        # 1. 检查起点 L+2 是否耗尽
        if L + 2 <= max_idx:
            if self._check_substrate_depletion(absorbances[L + 2], limit, direction):
                # L+2 点已耗尽，说明反应开始就没底物了
                return L, L, "No Linear Interval"

        # 2. 检查终点 M 是否耗尽，如果耗尽则向前查找 M'
        # 逻辑：从 M 开始往回找，直到找到一个没耗尽的点
        final_M = curr_M

        # 如果 limit 存在，且当前点耗尽，则缩短区间
        while final_M > L:
            if self._check_substrate_depletion(absorbances[final_M], limit, direction):
                final_M -= 1
            else:
                break  # 找到未耗尽点 M'

        return L, final_M, "OK"

    def _linearity_check(self, absorbances, L, M, total_slope, cycle_time, linearity_limit=0.2):
        """
        【4.5.5 线性度判断】
        公式: |Af - Ab| / |Atotal| < Limit
        """
        points_count = M - L + 1

        # 阈值检查：如果总反应率太小 (小于 0.006 Abs/min)，不做检查
        if abs(total_slope) <= self.MIN_slope_CHECK_THRESHOLD:
            return "Pass (Slope Low)"

        indices = list(range(L, M + 1))

        # 根据点数 N 决定分段方式
        if points_count > 8:
            # 前6个点，后6个点
            idx_f = indices[:6]
            idx_b = indices[-6:]
        elif 4 <= points_count <= 8:
            # 前3个点，后3个点
            idx_f = indices[:3]
            idx_b = indices[-3:]
        else:
            # N < 4，不做线性检查
            return "Skipped (N<4)"

        # 计算前后段斜率
        slope_f = self._calculate_slope_least_squares(absorbances, idx_f, cycle_time)
        slope_b = self._calculate_slope_least_squares(absorbances, idx_b, cycle_time)

        # 计算分子分母
        numerator = abs(slope_f - slope_b)
        denominator = abs(total_slope)

        # 分子阈值检查 (文档提到的 |Af - Ab| <= 60单位)
        if numerator <= self.MIN_slope_CHECK_THRESHOLD:
            return "Pass (Diff Low)"

        linearity = numerator / denominator

        if linearity > linearity_limit:
            return f"Fail (LIN) val={linearity:.2f}"

        return "Pass"

    def _enzyme_expansion(self, absorbances, delay_start_idx, L, cycle_time, limit, direction):
        """
        【4.5.6 酶线性范围扩展】
        当正常计算失败(点数不够)时触发。
        在延迟时间段 (delay_start ~ L) 寻找最大的两点间速率。
        """
        # 1. 确定延迟时间内的搜索范围 t1 ~ tL'
        # 也就是在延迟期内，也要排除掉底物耗尽的点
        search_end = L
        while search_end > delay_start_idx:
            if self._check_substrate_depletion(absorbances[search_end], limit, direction):
                search_end -= 1
            else:
                break

        valid_points = search_end - delay_start_idx + 1
        if valid_points < 2:
            return 0.0, "ENC (Expansion Calc Fail)"  # 无计算区间报错

        # 2. 遍历该区间，计算每两个相邻点的速率，取最大值 (Abs/min)
        # 公式: Delta A = 60 * (Ai+1 - Ai) / (ti+1 - ti)
        max_slope = 0.0
        found = False

        for i in range(delay_start_idx, search_end):
            val_curr = absorbances[i]
            val_next = absorbances[i + 1]

            # 计算两点斜率
            slope = 60 * (val_next - val_curr) / cycle_time

            # 比较绝对值取最大 (因为可能是负反应)
            if not found or abs(slope) > abs(max_slope):
                max_slope = slope
                found = True

        return max_slope, "EXP"

    def calculate_response(self, absorbances, N, P, L, M, K, cycle_time=1.0,
                           depletion_limit=None, direction=1, lin_limit=0.2):
        """
        主计算函数
        :param absorbances: 吸光度数组
        :param N, P: 空白区间索引
        :param L, M: 反应区间索引
        :param K: 空白扣除系数
        :param cycle_time: 周期时间(s)
        :param depletion_limit: 底物耗尽限值 (Abs)
        :param direction: 反应方向 1(增) / -1(减)
        :param lin_limit: 线性度限值 (默认0.2)
        """

        # --- 1. 确定线性范围 (4.5.3) ---
        real_L, real_M, range_status = self._determine_linear_range(absorbances, L, M, depletion_limit, direction)

        if range_status == "No Linear Interval":
            return 0.0, "Error: Substrate Depleted Early"

        points_count = real_M - real_L + 1

        # --- 2. 判断是否需要酶线性扩展 (4.5.6) ---
        # 若有效点数 < 2 (N=0 or 1)，启动扩展
        if points_count < 2:
            # 假设延迟时间是从索引 0 到 L (或者根据文档 t1 通常指反应开始)
            # 这里假设 t1 是 0，或者是 N(如果存在)
            delay_start = 0
            if P > 0: delay_start = P  # 如果有空白期，延迟期从空白期后算起

            exp_rate, exp_flag = self._enzyme_expansion(absorbances, delay_start, L, cycle_time, depletion_limit,
                                                        direction)
            return exp_rate, exp_flag

        # --- 3. 计算主反应速率 Delta A_LM (4.5.4) ---
        # 构造索引列表
        reaction_indices = list(range(real_L, real_M + 1))
        slope_reaction = self._calculate_slope_least_squares(absorbances, reaction_indices, cycle_time)

        # --- 4. 计算空白速率 Delta A_NP ---
        slope_blank = 0.0
        if N == 0 and P == 0:
            slope_blank = 0.0
        elif P > N:
            blank_indices = list(range(N, P + 1))
            slope_blank = self._calculate_slope_least_squares(absorbances, blank_indices, cycle_time)

        # --- 5. 计算最终反应度 R ---
        # R = Delta A_LM - K * Delta A_NP
        R = slope_reaction - (K * slope_blank)

        # --- 6. 线性度判断 (4.5.5) ---
        # 只有当点数足够 (N>=3, 实际上Check函数里分段需要N>=4) 且未触发扩展时进行
        lin_status = "Skipped"
        if points_count >= 4:
            lin_status = self._linearity_check(absorbances, real_L, real_M, slope_reaction, cycle_time, lin_limit)

            # 如果线性检查失败，系统应当打上 "LIN" 标记
            if "Fail" in lin_status:
                return R, f"Result: {R:.4f} (Flag: LIN)"

        # 若点数为 2，根据文档 4.5.3，给出报警但计算结果
        if points_count == 2:
            return R, f"Result: {R:.4f} (Flag: Points<3 Warning)"

        return R, f"Result: {R:.4f} (Linearity: {lin_status})"


class RateMethodCalculator:
    """
    【模块：速率法 / 动力学法计算】
    (Kinetic / Rate Method)
    通常用于酶类项目 (ALT, AST 等)。
    核心算法：最小二乘法线性拟合 (Least Squares Linear Regression)
    """

    def _linear_regression(self, x, y):
        """
        辅助函数：手写最小二乘法计算斜率
        y = kx + b, 返回 k
        """
        n = len(x)
        if n < 2: return 0.0

        sum_x = sum(x)
        sum_y = sum(y)
        sum_xy = sum(xi * yi for xi, yi in zip(x, y))
        sum_xx = sum(xi * xi for xi in x)

        numerator = n * sum_xy - sum_x * sum_y
        denominator = n * sum_xx - sum_x ** 2

        if denominator == 0:
            return 0.0

        slope = numerator / denominator
        return slope

    def calculate_response(self, absorbances, start_idx, end_idx, cycle_time=1.0):
        """
        计算速率法的反应度

        :param absorbances: 吸光度列表
        :param start_idx: 计算开始点 (L)
        :param end_idx: 计算结束点 (M)
        :param cycle_time: 测光周期(秒)
        :return: 反应度 R (Delta Abs / min)
        """

        # 1. 截取需要计算的数据段
        # 注意：生化仪通常会忽略前面不稳定的点，只取线性段
        if start_idx >= len(absorbances) or end_idx >= len(absorbances):
            print("[错误] 索引超出范围")
            return 0.0

        # 切片数据 (包含 end_idx)
        y_data = absorbances[start_idx: end_idx + 1]
        n_points = len(y_data)

        if n_points < 3:
            print("[警告] 点数过少，建议使用至少3个点进行拟合")
            # 如果点太少，退化为两点法计算
            if n_points == 2:
                slope = (y_data[-1] - y_data[0]) / cycle_time
                return slope * 60
            return 0.0

        # 2. 构建 X 轴 (时间)
        # 为了数值稳定性，通常从0开始构建相对时间
        # x = [0, cycle_time, 2*cycle_time, ...]
        x_data = [i * cycle_time for i in range(n_points)]

        # 3. 使用最小二乘法计算斜率 (Slope)
        # Slope 的单位是 Abs/sec
        slope = self._linear_regression(x_data, y_data)

        # 4. 转换为每分钟变化率
        R = slope * 60.0

        return R, slope, n_points

    def linearity_check(self, absorbances, start_idx, end_idx):
        """
        (进阶功能) 线性度检查
        有些生化仪要求检查反应曲线是否真的为直线。
        如果曲线弯曲（如底物耗尽），需要报错或缩短计算区间。
        """
        # 简单实现：将数据分为前后两半，分别算斜率，看差异是否过大
        mid = (start_idx + end_idx) // 2
        if mid <= start_idx: return True  # 点太少，无法检查

        # 计算前半段斜率
        # 注意：这里没传 cycle_time，只比较相对斜率即可
        calc = RateMethodCalculator()  # 实例化自身或调用静态方法
        r1, _, _ = self.calculate_response(absorbances, start_idx, mid)
        r2, _, _ = self.calculate_response(absorbances, mid, end_idx)

        # 如果前后斜率偏差超过 10% (0.1)，认为非线性
        if abs(r1) > 0.001:  # 防止分母为0
            diff = abs(r1 - r2) / abs(r1)
            if diff > 0.1:
                print(f"[非线性警告] 前段速率:{r1:.4f}, 后段速率:{r2:.4f}, 偏差:{diff * 100:.1f}%")
                return False
        return True


class EndpointCalculator:
    """
    【模块：终点法反应度计算】
    基于文档 4.3.5 和 4.3.6 实现。
    核心公式: R = Ai - k * Ab
    """

    def _calculate_mean_abs(self, absorbances, start_idx, end_idx):
        """
        辅助函数：计算指定区间 [start, end] 的吸光度平均值
        注意：文档中的 N<=P 和 L<=M 暗示这是一个闭区间
        """
        if not absorbances:
            return 0.0

        # 边界检查
        max_len = len(absorbances)
        if start_idx >= max_len:
            return 0.0

        # 修正结束索引，防止越界
        real_end = min(end_idx, max_len - 1)

        if start_idx > real_end:
            return 0.0

        # 切片并计算平均值
        # Python切片是左闭右开，所以 end_idx 需要 +1
        segment = absorbances[start_idx: real_end + 1]
        return sum(segment) / len(segment)

    def calculate_raw_response(self, absorbances, N, P, L, M, K):
        """
        计算基础反应度 R
        对应文档 4.3.5

        :param absorbances: 反应曲线数据 (list)
        :param N: 空白时间起始点 (索引)
        :param P: 空白时间结束点 (索引)
        :param L: 反应时间起始点 (索引)
        :param M: 反应时间结束点 (索引)
        :param K: 体积校正因子 (k)
        :return: 计算出的 R 值
        """
        # 1. 计算反应时间吸光度 Ai (L ~ M)
        Ai = self._calculate_mean_abs(absorbances, L, M)

        # 2. 计算空白时间吸光度 Ab (N ~ P)
        # 表 4.1 中提到 "不扣除空白" 时 N=P=0，且 K 通常为 0，
        # 但为了逻辑严谨，若 K=0，其实 Ab 算不算都行，这里照常计算。
        if N == 0 and P == 0 and K == 0:
            Ab = 0.0
        else:
            Ab = self._calculate_mean_abs(absorbances, N, P)

        # 3. 应用公式 R = Ai - k * Ab
        R = Ai - (K * Ab)

        return R, Ai, Ab

    def calculate_corrected_response(self,
                                     sample_data,
                                     sample_blank_data,
                                     N, P, L, M, K):
        """
        计算样本空白校正后的反应度 R'
        对应文档 4.3.6
        公式: R' = R - Rsb

        :param sample_data: 样本反应曲线
        :param sample_blank_data: 样本空白反应曲线 (如无则传 None)
        :param N, P, L, M, K: 同上，计算参数
        """
        # 1. 计算样本本身的 R
        R_sample, _, _ = self.calculate_raw_response(sample_data, N, P, L, M, K)

        # 2. 如果没有样本空白数据，直接返回 R
        if not sample_blank_data:
            return R_sample

        # 3. 计算样本空白的 Rsb
        # 文档：样本空白反应度也以上述 R 公式计算
        R_sb, _, _ = self.calculate_raw_response(sample_blank_data, N, P, L, M, K)

        # 4. 计算校正后结果 R' = R - Rsb
        R_prime = R_sample - R_sb

        print(f"[计算细节] R样本={R_sample:.4f}, R样空={R_sb:.4f}, R'={R_prime:.4f}")
        return R_prime


class LinearCalibrator:
    """
    【模块：线性校准 (Linear Calibration)】
    基于文档 4.6.1 实现。
    核心公式: C = K * (R - R0)
    关键处理: R 和 R0 需要除以 10000 (还原吸光度单位)
    """

    def __init__(self):
        self.K = 1.0
        self.R0 = 0.0
        self.is_calibrated = False
        self.SCALE_FACTOR = 10000.0  # 文档规定的除数

    def _scale_R(self, raw_R):
        """
        内部辅助：处理除以 10000 的逻辑
        """
        return raw_R / self.SCALE_FACTOR

    def fit_single_point(self, user_K, user_R0_raw=0.0):
        """
        1. 单点线性校准 (K因数法)
        :param user_K: 用户输入的 K 因数
        :param user_R0_raw: 试剂空白反应度 (未除以10000前的原始值)，若无则为0
        """
        self.K = user_K
        # R0 也需要除以 10000
        self.R0 = self._scale_R(user_R0_raw)
        self.is_calibrated = True
        print(f"[校准成功] 单点法: K={self.K:.4f}, R0={self.R0:.4f}")

    def fit_two_point(self, concs, resps_raw):
        """
        2. 两点线性校准
        要求提供 2 个校准品。
        """
        if len(concs) != 2 or len(resps_raw) != 2:
            print("[错误] 两点法必须传入 2 个点的浓度和反应度")
            return

        # 1. 数据准备与缩放
        C1, C2 = concs[0], concs[1]
        # 将原始反应度除以 10000
        R1 = self._scale_R(resps_raw[0])
        R2 = self._scale_R(resps_raw[1])

        # 2. 防止分母为0
        if abs(R2 - R1) < 1e-9:
            print("[错误] 两个校准品的反应度相同，无法计算斜率")
            return

        # 3. 计算 K (斜率)
        # 公式: K = (C2 - C1) / (R2 - R1)
        self.K = (C2 - C1) / (R2 - R1)

        # 4. 计算 R0 (截距相关项)
        # 公式: R0 = R1 - C1 / K
        if self.K == 0:
            self.R0 = 0
        else:
            self.R0 = R1 - (C1 / self.K)

        self.is_calibrated = True
        print(f"[校准成功] 两点法: K={self.K:.4f}, R0={self.R0:.4f}")

    def fit_multi_point(self, concs, resps_raw):
        """
        3. 多点线性校准 (最小二乘法)
        要求提供 n (>=3) 个校准品。
        """
        n = len(concs)
        if n < 3:
            print("[错误] 多点法至少需要 3 个校准品")
            return
        if len(concs) != len(resps_raw):
            print("[错误] 浓度与反应度数量不匹配")
            return

        # 1. 数据准备 (应用除以 10000 规则)
        C = concs
        R = [self._scale_R(r) for r in resps_raw]

        # 2. 计算求和项
        sum_C = sum(C)
        sum_R = sum(R)

        sum_CR = sum(c * r for c, r in zip(C, R))  # sum(Ci * Ri)
        sum_R2 = sum(r * r for r in R)  # sum(Ri^2)

        # 3. 计算 K
        # 分子: sum(Ci*Ri) - (sum(Ci)*sum(Ri))/n
        numerator = sum_CR - (sum_C * sum_R) / n

        # 分母: sum(Ri^2) - (sum(Ri)^2)/n
        denominator = sum_R2 - (sum_R ** 2) / n

        if abs(denominator) < 1e-9:
            print("[错误] 反应度数据方差为0，无法拟合")
            return

        self.K = numerator / denominator

        # 4. 计算 R0
        # 公式: R0 = (sum(Ri)/n) - (sum(Ci)/n)/K
        # 其实就是: R0 = mean(R) - mean(C)/K
        # 这意味着标准方程 C = K*R - K*R0  => mean(C) = K*mean(R) - K*R0 => K*R0 = K*mean(R) - mean(C)
        if self.K == 0:
            self.R0 = 0
        else:
            mean_R = sum_R / n
            mean_C = sum_C / n
            self.R0 = mean_R - (mean_C / self.K)

        self.is_calibrated = True
        print(f"[校准成功] 多点法({n}点): K={self.K:.4f}, R0={self.R0:.4f}")

    def calculate_concentration(self, sample_resp_raw):
        """
        计算浓度
        公式: C = K * (R - R0)
        注意: 输入的 sample_resp_raw 也会在此处被除以 10000
        """
        if not self.is_calibrated:
            print("[警告] 尚未校准")
            return 0.0

        # 1. 对样本反应度进行缩放
        R = self._scale_R(sample_resp_raw)

        # 2. 代入公式
        C = self.K * (R - self.R0)

        # 浓度通常不为负，除非做偏移修正，这里保留原始计算结果
        return C


# ====================================================================
#
#                生化分析仪核心计算方法演示
#
# ====================================================================

if __name__ == "__main__":
    # ========================================================
    # 模块一: 固定时间法 (Fixed-Time Method)
    # ========================================================
    print("\n" + "="*50)
    print("模块一: 固定时间法 (Fixed-Time Method)")
    print("="*50)
    print("适用于反应速率在监测时间内恒定的项目，如肌酐(苦味酸法)。\n")

    # --- 数据准备 ---
    CYCLE_TIME_SEC = 18.0
    calc = FixedTimeCalculator()
    # 模拟数据：苦味酸法测肌酐 (Jaffe Method)
    raw_curve = []
    # 0-16点：试剂1 + 样本 (空白期)
    for i in range(17):
        raw_curve.append(0.10 + i * 0.001)
    # 17-33点：加入试剂2 (主反应期)
    for i in range(20):
        raw_curve.append(0.12 + i * 0.02)

    print(f"-> 模拟数据: 肌酐反应曲线 (Jaffe Method)")
    print(f"  - 数据点总数: {len(raw_curve)}")
    print(f"  - 测光周期: {CYCLE_TIME_SEC} 秒\n")

    # --- 场景 1 ---
    print("--- 场景 1: 双试剂，双区间扣除 (带体积修正) ---")
    print("   目的: 扣除主反应前由试剂1引起的非特异性反应速率。")
    N, P = 5, 16
    L, M = 17, 30
    K_factor = 0.8
    R, s_react, s_blank = calc.calculate_response(raw_curve, N, P, L, M, K_factor, CYCLE_TIME_SEC)

    print(f"   输入参数:")
    print(f"     - 空白监测窗 (N-P): {N}-{P}")
    print(f"     - 反应监测窗 (L-M): {L}-{M}")
    print(f"     - 体积修正因子 (K): {K_factor}")
    print(f"   计算过程:")
    print(f"     - 反应区间斜率: {s_react:.6f} Abs/sec")
    print(f"     - 空白区间斜率: {s_blank:.6f} Abs/sec")
    print(f"     - 计算公式: R = 60 * (反应斜率 - K * 空白斜率)")
    print(f"     -           = 60 * ({s_react:.6f} - {K_factor} * {s_blank:.6f})")
    print(f"   最终结果:")
    print(f"     - 最终反应度 R: {R:.4f} ΔAbs/min\n")

    # --- 场景 2 ---
    print("--- 场景 2: 不扣除空白 (K=0) ---")
    print("   目的: 仅计算主反应区间的速率，适用于无明显试剂空白的项目。")
    N, P = 0, 0
    L, M = 17, 30
    K_factor = 0.0
    R, s_react, s_blank = calc.calculate_response(raw_curve, N, P, L, M, K_factor, CYCLE_TIME_SEC)

    print(f"   输入参数:")
    print(f"     - 空白监测窗 (N-P): {N}-{P} (不使用)")
    print(f"     - 反应监测窗 (L-M): {L}-{M}")
    print(f"     - 体积修正因子 (K): {K_factor}")
    print(f"   计算过程:")
    print(f"     - 反应区间斜率: {s_react:.6f} Abs/sec")
    print(f"     - 计算公式: R = 60 * 反应斜率")
    print(f"   最终结果:")
    print(f"     - 最终反应度 R: {R:.4f} ΔAbs/min\n")


    # ========================================================
    # 模块二: 终点法 (Endpoint Method)
    # ========================================================
    print("\n" + "="*50)
    print("模块二: 终点法 (Endpoint Method)")
    print("="*50)
    print("适用于反应达到或接近终点时进行测量的项目，如总蛋白、白蛋白等。\n")

    # --- 数据准备 ---
    calc = EndpointCalculator()
    raw_curve = [0.1] * 10 + [0.12] * 30 + [0.5 + 0.01 * i for i in range(30)]
    print(f"-> 模拟数据: 一个典型的终点法反应曲线")
    print(f"  - 数据点总数: {len(raw_curve)}\n")


    # --- 场景 1 ---
    print("--- 场景 1: 双试剂，空白时间在反应启动前 (带体积修正) ---")
    N, P = 5, 15
    L, M = 50, 60
    K_factor = 0.8
    R, Ai, Ab = calc.calculate_raw_response(raw_curve, N, P, L, M, K_factor)

    print(f"   输入参数:")
    print(f"     - 空白窗 (N-P): {N}-{P}")
    print(f"     - 反应窗 (L-M): {L}-{M}")
    print(f"     - 体积修正因子 (K): {K_factor}")
    print(f"   计算过程:")
    print(f"     - 反应终点均值 (Ai): {Ai:.4f}")
    print(f"     - 空白读数均值 (Ab): {Ab:.4f}")
    print(f"     - 计算公式: R = Ai - K * Ab")
    print(f"   最终结果:")
    print(f"     - 最终反应度 R = {Ai:.4f} - {K_factor} * {Ab:.4f} = {R:.4f}\n")


    # --- 场景 2 ---
    print("--- 场景 2: 单试剂，扣除初期读数 (K=1) ---")
    N, P = 5, 6
    L, M = 60, 65
    K_factor = 1.0
    R, Ai, Ab = calc.calculate_raw_response(raw_curve, N, P, L, M, K_factor)

    print(f"   输入参数:")
    print(f"     - 空白窗 (N-P): {N}-{P} (反应初期)")
    print(f"     - 反应窗 (L-M): {L}-{M} (反应终点)")
    print(f"     - K 因子: {K_factor} (表示直接相减)")
    print(f"   最终结果:")
    print(f"     - 最终反应度 R = {Ai:.4f} - {K_factor:.1f} * {Ab:.4f} = {R:.4f}\n")

    # --- 场景 3 ---
    print("--- 场景 3: 不扣除空白 (K=0) ---")
    N, P = 0, 0
    L, M = 60, 65
    K_factor = 0.0
    R, Ai, Ab = calc.calculate_raw_response(raw_curve, N, P, L, M, K_factor)

    print(f"   输入参数:")
    print(f"     - 空白窗 (N-P): {N}-{P} (不使用)")
    print(f"     - 反应窗 (L-M): {L}-{M}")
    print(f"     - K 因子: {K_factor} (不扣除空白)")
    print(f"   最终结果:")
    print(f"     - 最终反应度 R = {Ai:.4f}\n")

    # --- 场景 4 ---
    print("--- 场景 4: 样本空白校正 ---")
    print("   目的: 扣除由样本自身颜色、浊度等引起的背景吸收。")
    sample_curve = [x + 0.5 for x in raw_curve]
    sample_blank_curve = [0.52] * 70
    N, P, L, M, K = 5, 15, 50, 60, 0.8

    final_R = calc.calculate_corrected_response(sample_curve, sample_blank_curve, N, P, L, M, K)
    print(f"   说明: 使用场景1的参数，对一个背景增高(如脂血)的样本进行校正。")
    print(f"   最终结果:")
    print(f"     - 校正后反应度 R': {final_R:.4f}\n")


    # ========================================================
    # 模块三: 动力学法 (Kinetic Method)
    # ========================================================
    print("\n" + "="*50)
    print("模块三: 动力学法 (Kinetic Method)")
    print("="*50)
    print("主要用于酶活性检测，通过监测反应速率来计算。包含复杂的线性检查和范围调整逻辑。\n")

    # --- 数据准备 ---
    calc = KineticCalculator()
    CYCLE_TIME = 10.0
    raw_data = [2.5, 2.2, 1.8, 1.4, 1.0, 0.6, 0.4, 0.3, 0.28, 0.28, 0.28, 0.28, 0.28]
    N, P, K = 0, 0, 0
    LIMIT = 0.5
    DIR = -1
    print(f"-> 模拟数据: 高活性酶反应，导致底物耗尽")
    print(f"  - 反应方向: 下降 (DIR={DIR})")
    print(f"  - 耗尽限值: {LIMIT} Abs")
    print(f"  - 测光周期: {CYCLE_TIME} 秒\n")

    # --- 测试 1 ---
    print("--- 测试 1: 底物耗尽，触发线性范围缩减 ---")
    L_test, M_test = 3, 10
    R, msg = calc.calculate_response(raw_data, N, P, L_test, M_test, K, CYCLE_TIME, LIMIT, DIR)
    print(f"   输入参数:")
    print(f"     - 预设区间 (L-M): [{L_test}~{M_test}]")
    print(f"   运行信息: {msg}")
    print(f"   说明: 算法检测到在 M={M_test} 点时反应已耗尽，自动缩短计算区间以获得有效速率。\n")

    # --- 测试 2 ---
    print("--- 测试 2: 起点耗尽，触发酶线性扩展 (EXP) ---")
    L, M = 6, 12
    R, msg = calc.calculate_response(raw_data, N, P, L, M, K, CYCLE_TIME, LIMIT, DIR)
    print(f"   输入参数:")
    print(f"     - 预设区间 (L-M): [{L}~{M}]")
    print(f"   运行信息: {msg}")
    print(f"   说明: 算法检测到在 L={L} 点时反应就已耗尽，预设主区间完全无效。")
    print(f"           自动回溯到延迟期，寻找最快的反应速率作为结果，并标记为 'EXP'。\n")

    # --- 测试 3 ---
    print("--- 测试 3: 正常线性反应 (通过线性度检查) ---")
    linear_data = [1.0 - 0.01 * i for i in range(20)]
    L_lin, M_lin = 5, 15
    R, msg = calc.calculate_response(linear_data, N, P, L_lin, M_lin, K, CYCLE_TIME, None, DIR)
    print(f"   输入参数:")
    print(f"     - 预设区间 (L-M): [{L_lin}~{M_lin}]")
    print(f"   运行信息: {msg}")
    print(f"   说明: 在理想的线性数据上，算法通过线性度检查，正常计算结果。\n")

    # --- 测试 4 ---
    print("--- 测试 4: 非线性反应 (触发 LIN 报警) ---")
    curved_data = [1.0 - 0.02 * i for i in range(10)] + [0.8 - 0.005 * i for i in range(10)]
    L_cur, M_cur = 0, 18
    R, msg = calc.calculate_response(curved_data, N, P, L_cur, M_cur, K, CYCLE_TIME, None, DIR)
    print(f"   输入参数:")
    print(f"     - 预设区间 (L-M): [{L_cur}~{M_cur}]")
    print(f"   运行信息: {msg}")
    print(f"   说明: 数据曲线发生弯曲 (速率变化)，算法的线性度检查失败，并标记 'LIN' 报警。\n")


    # ========================================================
    # 模块四: 速率法 (Rate Method)
    # ========================================================
    print("\n" + "="*50)
    print("模块四: 速率法 (Rate Method)")
    print("="*50)
    print("是动力学法的一种简化形式，通常指用最小二乘法对指定区间的数据进行线性拟合来求速率。\n")

    # --- 数据准备 ---
    calc = RateMethodCalculator()
    import random
    raw_curve = [(1.5 + (-0.002 * i * 18) + random.uniform(-0.0005, 0.0005)) for i in range(50)]
    L, M = 10, 30
    CYCLE_TIME = 18.0
    true_slope_per_sec = -0.002

    print(f"-> 模拟数据: ALT (谷丙转氨酶) 反应曲线 (带随机噪声)")
    print(f"  - 测光周期: {CYCLE_TIME} 秒")
    print(f"  - 理论斜率: {true_slope_per_sec} Abs/sec\n")

    # --- 计算 ---
    print("--- 速率法计算演示 ---")
    R, slope, pts = calc.calculate_response(raw_curve, L, M, CYCLE_TIME)
    is_linear = calc.linearity_check(raw_curve, L, M)

    print(f"   输入参数:")
    print(f"     - 计算区间 (L-M): {L}-{M}")
    print(f"   计算结果:")
    print(f"     - 参与拟合点数: {pts}")
    print(f"     - 最小二乘法拟合斜率: {slope:.6f} Abs/sec")
    print(f"     - 最终反应度 R: {R:.4f} ΔAbs/min")
    print(f"   对比:")
    print(f"     - 理论真值: {true_slope_per_sec * 60:.4f} ΔAbs/min")
    print(f"   检查项:")
    print(f"     - 线性检查结果: {'通过' if is_linear else '失败'}\n")


    # ========================================================
    # 模块五: 线性校准 (Linear Calibration)
    # ========================================================
    print("\n" + "="*50)
    print("模块五: 线性校准 (Linear Calibration)")
    print("="*50)
    print("通过标准品建立反应度(R)与浓度(C)之间的线性关系 C = a*R + b。\n")

    calibrator = LinearCalibrator()

    # --- 测试 1 ---
    print("--- 测试 1: 单点线性校准 (K因数法) ---")
    calibrator.fit_single_point(user_K=100.0, user_R0_raw=500)
    sample_val = 2000
    conc = calibrator.calculate_concentration(sample_val)
    print(f"   校准参数: K = 100.0, R0 = 500 (原始信号值)")
    print(f"   样本计算: 原始反应度 {sample_val} -> 计算浓度: {conc:.4f}\n")

    # --- 测试 2 ---
    print("--- 测试 2: 两点线性校准 ---")
    raw_resps = [1000, 5000]
    std_concs = [10.0, 50.0]
    calibrator.fit_two_point(std_concs, raw_resps)
    sample_val = 2000
    conc = calibrator.calculate_concentration(sample_val)
    print(f"   校准点: C1={std_concs[0]} @ R1={raw_resps[0]}, C2={std_concs[1]} @ R2={raw_resps[1]}")
    print(f"   样本计算: 原始反应度 {sample_val} -> 计算浓度: {conc:.4f}\n")

    # --- 测试 3 ---
    print("--- 测试 3: 多点线性校准 (最小二乘法) ---")
    noisy_resps = [1000, 3000, 5000, 7000]
    noisy_concs = [10.0, 30.0, 50.0, 68.0]
    calibrator.fit_multi_point(noisy_concs, noisy_resps)
    test_val = 4000
    conc = calibrator.calculate_concentration(test_val)
    print(f"   校准点 (带噪声): C={noisy_concs} @ R={noisy_resps}")
    print(f"   样本计算: 原始反应度 {test_val} -> 计算浓度: {conc:.4f}\n")


    # ========================================================
    # 模块六: 非线性校准 (Non-Linear Calibration)
    # ========================================================
    print("\n" + "="*50)
    print("模块六: 非线性校准 (Non-Linear Calibration)")
    print("="*50)
    print("用于免疫等反应，其反应度(R)与浓度(C)不成线性关系，需要用复杂曲线拟合。\n")
    calib = NonLinearCalibrator()

    # --- 测试 1 ---
    print("--- 测试 1: Polynomial 5P 模型 ---")
    concs = [0, 10, 50, 100, 200]
    resps = [500, 1200, 2500, 3200, 3800]
    calib.fit_polynomial_5p(concs, resps)
    sample_R = 2000
    conc = calib.calculate_concentration(sample_R)
    print(f"   校准点: C={concs} @ R={resps}")
    print(f"   样本计算: R={sample_R} -> Conc={conc:.4f} (预期在 10-50 之间)\n")

    # --- 测试 2 ---
    print("--- 测试 2: Parabola (二次多项式) 模型 ---")
    c_para = [0, 50, 100]
    r_para = [100, 5000, 9000]
    calib.fit_parabola(c_para, r_para)
    conc_p = calib.calculate_concentration(3000)
    print(f"   校准点: C={c_para} @ R={r_para}")
    print(f"   样本计算: R=3000 -> Conc={conc_p:.4f}\n")

    # --- 测试 3 ---
    print("--- 测试 3: Logistic 4P 模型 (直接使用参数) ---")
    calib.set_params_logistic_4p(R0=100, K=5000, a=2.0, b=-1.0)
    conc_l = calib.calculate_concentration(2600)
    print(f"   模型参数: R0=100, K=5000, a=2.0, b=-1.0")
    print(f"   样本计算: R=2600 -> Conc={conc_l:.4f}\n")
