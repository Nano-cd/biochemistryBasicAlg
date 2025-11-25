#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <iomanip>
#include <utility> // For std::pair
#include <sstream> // <--- 添加这一行

// 类型别名，使代码更清晰
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

// ==========================================
// 基础线性代数工具 (用于替代 numpy)
// ==========================================
/**
 * @brief 基础线性代数工具类
 * @details 提供矩阵转置、乘法、求逆和最小二乘法求解等静态方法。
 */
class MatrixUtils {
public:
	static Matrix transpose(const Matrix& matrix) {
		if (matrix.empty()) return {};
		size_t m = matrix.size();
		size_t n = matrix[0].size();
		Matrix result(n, Vector(m));
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < n; ++j) {
				result[j][i] = matrix[i][j];
			}
		}
		return result;
	}

	static Matrix multiply(const Matrix& A, const Matrix& B) {
		if (A.empty() || B.empty()) return {};
		size_t m = A.size();
		size_t n = A[0].size();
		size_t p = B[0].size();
		if (n != B.size()) {
			throw std::invalid_argument("矩阵维度不匹配，无法相乘");
		}
		Matrix C(m, Vector(p, 0.0));
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < p; ++j) {
				for (size_t k = 0; k < n; ++k) {
					C[i][j] += A[i][k] * B[k][j];
				}
			}
		}
		return C;
	}

	static Matrix inverse(Matrix A) {
		size_t n = A.size();
		if (n == 0 || n != A[0].size()) {
			throw std::invalid_argument("输入必须为方阵");
		}

		Matrix M(n, Vector(2 * n));
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				M[i][j] = A[i][j];
			}
			M[i][i + n] = 1.0;
		}

		for (size_t i = 0; i < n; ++i) {
			double pivot = M[i][i];
			if (std::abs(pivot) < 1e-12) {
				throw std::runtime_error("矩阵不可逆");
			}
			for (size_t j = 0; j < 2 * n; ++j) {
				M[i][j] /= pivot;
			}
			for (size_t k = 0; k < n; ++k) {
				if (k != i) {
					double factor = M[k][i];
					for (size_t j = 0; j < 2 * n; ++j) {
						M[k][j] -= factor * M[i][j];
					}
				}
			}
		}

		Matrix inv(n, Vector(n));
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j) {
				inv[i][j] = M[i][j + n];
			}
		}
		return inv;
	}

	static Vector solve_least_squares(const Matrix& X, const Matrix& Y) {
		Matrix XT = transpose(X);
		Matrix XTX = multiply(XT, X);
		try {
			Matrix XTX_inv = inverse(XTX);
			Matrix XTY = multiply(XT, Y);
			Matrix Beta_mat = multiply(XTX_inv, XTY);

			Vector Beta;
			for (const auto& row : Beta_mat) {
				Beta.push_back(row[0]);
			}
			return Beta;
		}
		catch (const std::runtime_error&) {
			return {}; // 返回空向量表示失败
		}
	}
};

// ==========================================
// 非线性校准核心类
// ==========================================
/**
 * @brief 非线性校准核心类
 * @details 对应文档 4.6.2 非线性校准。
 *          实现了多种非线性模型，如 Polynomial 5P, Parabola, Logistic 4P, Spline 等。
 */
class NonLinearCalibrator {
private:
	std::string method;
	Vector params;
	std::vector<std::pair<double, double>> spline_data;

public:
	void fit_polynomial_5p(const Vector& concs, const Vector& resps) {
		if (concs.size() < 5) {
			std::cerr << "[Error] Polynomial 5P 需要至少 5 个校准品" << std::endl;
			return;
		}

		double R0 = resps[0];
		Matrix X_mat;
		Matrix Y_mat;
		int valid_points = 0;

		for (size_t i = 0; i < concs.size(); ++i) {
			if (concs[i] <= 0) continue;

			double x_val = (resps[i] - R0) / 100.0;
			double ln_c = std::log(concs[i]);
			X_mat.push_back({ 1.0, x_val, std::pow(x_val, 2), std::pow(x_val, 3) });
			Y_mat.push_back({ ln_c });
			valid_points++;
		}

		if (valid_points < 4) {
			std::cerr << "[Error] 有效点(C>0)不足以拟合3次方程" << std::endl;
			return;
		}

		Vector coeffs = MatrixUtils::solve_least_squares(X_mat, Y_mat);
		if (coeffs.empty()) {
			std::cerr << "[Error] 拟合失败" << std::endl;
			return;
		}

		this->method = "Polynomial_5P";
		this->params = { R0 };
		this->params.insert(this->params.end(), coeffs.begin(), coeffs.end());
		std::cout << "[拟合成功] Poly 5P: R0=" << R0 << ", a=" << std::fixed << std::setprecision(4) << coeffs[0] << ", b=" << coeffs[1] << "..." << std::endl;
	}

	double calc_polynomial_5p(double R) {
		double R0 = params[0], a = params[1], b = params[2], c = params[3], d = params[4];
		double x = (R - R0) / 100.0;
		double ln_c = a + b * x + c * std::pow(x, 2) + d * std::pow(x, 3);
		return std::exp(ln_c);
	}

	void fit_parabola(const Vector& concs, const Vector& resps) {
		if (concs.size() < 3) {
			std::cerr << "[Error] Parabola 需要至少 3 个校准品" << std::endl;
			return;
		}
		Matrix X_mat;
		Matrix Y_mat;
		for (size_t i = 0; i < concs.size(); ++i) {
			X_mat.push_back({ std::pow(concs[i], 2), concs[i], 1.0 });
			Y_mat.push_back({ resps[i] });
		}
		Vector coeffs = MatrixUtils::solve_least_squares(X_mat, Y_mat);
		this->method = "Parabola";
		this->params = coeffs;
		std::cout << "[拟合成功] Parabola: R=" << std::fixed << std::setprecision(4) << coeffs[0] << "C^2 + " << coeffs[1] << "C + " << coeffs[2] << std::endl;
	}

	double calc_parabola(double R) {
		double a = params[0], b = params[1], R0 = params[2];
		double delta = b * b - 4 * a * (R0 - R);
		if (delta < 0) return 0.0;
		double sqrt_delta = std::sqrt(delta);
		if (std::abs(a) < 1e-12) return (b != 0) ? (R - R0) / b : 0.0;

		double c1 = (-b + sqrt_delta) / (2 * a);
		double c2 = (-b - sqrt_delta) / (2 * a);

		if (c1 >= 0) return c1;
		if (c2 >= 0) return c2;
		return 0.0;
	}

	void set_params_logistic_4p(double R0, double K, double a, double b) {
		this->method = "Logistic_4P";
		this->params = { R0, K, a, b };
	}

	double calc_logistic_4p(double R) {
		double R0 = params[0], K = params[1], a = params[2], b = params[3];
		try {
			double val = K / (R - R0) - 1;
			if (val <= 0) return 0.0;
			double ln_E = std::log(val);
			double ln_C = -(ln_E + a) / b;
			return std::exp(ln_C);
		}
		catch (...) {
			return 0.0;
		}
	}

	void set_params_exponential_5p(double R0, double K, double a, double b, double c) {
		this->method = "Exponential_5P";
		this->params = { R0, K, a, b, c };
	}

	double calc_exponential_5p(double R) {
		double R0 = params[0], K = params[1], a = params[2], b = params[3], c_param = params[4];
		try {
			double target = std::log((R - R0) / K);
			double x = 1.0;
			for (int i = 0; i < 20; ++i) {
				double fx = a * x + b * x * x + c_param * x * x * x - target;
				double dfx = a + 2 * b * x + 3 * c_param * x * x;
				if (std::abs(dfx) < 1e-6) break;
				double x_new = x - fx / dfx;
				if (std::abs(x_new - x) < 1e-6) {
					x = x_new;
					break;
				}
				x = x_new;
			}
			return std::exp(x);
		}
		catch (...) {
			return 0.0;
		}
	}

	void fit_spline(const Vector& concs, const Vector& resps) {
		this->spline_data.clear();
		for (size_t i = 0; i < concs.size(); ++i) {
			this->spline_data.push_back({ resps[i], concs[i] });
		}
		std::sort(this->spline_data.begin(), this->spline_data.end());
		this->method = "Spline";
		std::cout << "[拟合成功] Spline: 存储了 " << this->spline_data.size() << " 个节点" << std::endl;
	}

	double calc_spline(double R) {
		if (spline_data.empty()) return 0.0;
		if (R <= spline_data.front().first) return spline_data.front().second;
		if (R >= spline_data.back().first) return spline_data.back().second;

		for (size_t i = 0; i < spline_data.size() - 1; ++i) {
			double r1 = spline_data[i].first, c1 = spline_data[i].second;
			double r2 = spline_data[i + 1].first, c2 = spline_data[i + 1].second;
			if (R >= r1 && R <= r2) {
				if (r2 - r1 == 0) return c1;
				return c1 + (R - r1) * (c2 - c1) / (r2 - r1);
			}
		}
		return 0.0;
	}

	double calculate_concentration(double R) {
		if (method == "Polynomial_5P") return calc_polynomial_5p(R);
		if (method == "Parabola") return calc_parabola(R);
		if (method == "Logistic_4P") return calc_logistic_4p(R);
		if (method == "Exponential_5P") return calc_exponential_5p(R);
		if (method == "Spline") return calc_spline(R);
		std::cerr << "未选择校准方法" << std::endl;
		return 0.0;
	}
};

/**
 * @brief 存储固定时间法计算结果的结构体
 */
struct FixedTimeResult {
    double R;               ///< 最终反应度 (ΔAbs/min)
    double slope_reaction;  ///< 主反应区间的斜率 (Abs/sec)
    double slope_blank;     ///< 空白区间的斜率 (Abs/sec)
};

/**
 * @brief 固定时间法反应度计算模块
 * @details 对应文档 4.4.1 和 4.4.2。又称一级动力学法或初速率法。
 *          适用于反应速率在监测时间内恒定的项目。
 */
class FixedTimeCalculator {
public:
	/**
 * @brief 计算固定时间法的反应度 R
 * @param[in] absorbances 吸光度数据列表
 * @param[in] N 空白区间起始点索引
 * @param[in] P 空白区间结束点索引
 * @param[in] L 反应区间起始点索引
 * @param[in] M 反应区间结束点索引
 * @param[in] K 体积校正因子
 * @param[in] cycle_time 测光周期(秒)
 * @return FixedTimeResult 包含最终反应度 R 和中间斜率值的结构体
 */
	FixedTimeResult calculate_response(const Vector& absorbances, int N, int P, int L, int M, double K, double cycle_time = 1.0) {
		size_t max_idx = absorbances.size() - 1;
		if (M > max_idx || L > max_idx) {
			std::cerr << "[错误] 反应区间索引超出数据范围" << std::endl;
			return { 0.0, 0.0, 0.0 };
		}

		double val_M = absorbances[M];
		double val_L = absorbances[L];
		double time_diff_reaction = (M - L) * cycle_time;

		if (time_diff_reaction == 0) {
			std::cerr << "[错误] 反应区间时间差为0 (L==M)" << std::endl;
			return { 0.0, 0.0, 0.0 };
		}
		double slope_reaction = (val_M - val_L) / time_diff_reaction;
		double slope_blank = 0.0;

		if (P > 0 && P > N && K != 0) {
			if (P > max_idx) {
				std::cerr << "[错误] 空白区间索引超出数据范围" << std::endl;
				return { 0.0, slope_reaction, 0.0 };
			}
			double val_P = absorbances[P];
			double val_N = absorbances[N];
			double time_diff_blank = (P - N) * cycle_time;
			if (time_diff_blank != 0) {
				slope_blank = (val_P - val_N) / time_diff_blank;
			}
		}
		double R = 60.0 * (slope_reaction - (K * slope_blank));
		return { R, slope_reaction, slope_blank };
	}
};

/**
 * @brief 动力学法反应度计算模块
 * @details 对应文档 4.5 全章节。主要用于酶活性检测。
 *          核心功能包括：线性范围动态判定、最小二乘法速率计算、线性度检查和酶线性范围扩展。
 */
class KineticCalculator {
private:
	/**
 * @brief 动力学法主计算函数
 * @param[in] absorbances 吸光度数据列表
 * @param[in] N 空白区间起始点索引
 * @param[in] P 空白区间结束点索引
 * @param[in] L 反应区间起始点索引
 * @param[in] M 反应区间结束点索引
 * @param[in] K 体积校正因子
 * @param[in] cycle_time 测光周期(秒)
 * @param[in] depletion_limit 底物耗尽限值 (吸光度单位)。如果吸光度超过此限值，则认为底物耗尽。
 * @param[in] direction 反应方向: 1 (增加) / -1 (减少)
 * @param[in] lin_limit 线性度检查的阈值，默认为 0.2 (20%)
 * @return std::pair<double, std::string> 返回计算出的反应度 R 和一个包含状态信息（如LIN, EXP报警）的字符串。
 */
	double MIN_slope_CHECK_THRESHOLD = 0.006;

	double _calculate_slope_least_squares(const Vector& absorbances, const std::vector<int>& indices, double cycle_time) {
		size_t n = indices.size();
		if (n < 2) return 0.0;

		Vector T, A;
		for (int idx : indices) {
			T.push_back(idx * cycle_time);
			A.push_back(absorbances[idx]);
		}

		double T_bar = std::accumulate(T.begin(), T.end(), 0.0) / n;
		double A_bar = std::accumulate(A.begin(), A.end(), 0.0) / n;

		double numerator = 0.0;
		double denominator = 0.0;
		for (size_t i = 0; i < n; ++i) {
			numerator += (T[i] - T_bar) * (A[i] - A_bar);
			denominator += (T[i] - T_bar) * (T[i] - T_bar);
		}
		return (denominator == 0) ? 0.0 : 60.0 * (numerator / denominator);
	}

	bool _check_substrate_depletion(double abs_val, double limit, int direction) {
		if (direction == 1 && abs_val > limit) return true;
		if (direction == -1 && abs_val < limit) return true;
		return false;
	}

	struct LinearRangeResult { int L, M; std::string status; };
	LinearRangeResult _determine_linear_range(const Vector& absorbances, int L, int M, double depletion_limit, int direction) {
		size_t max_idx = absorbances.size() - 1;
		int curr_M = std::min((size_t)M, max_idx);

		if (L + 2 <= max_idx) {
			if (_check_substrate_depletion(absorbances[L + 2], depletion_limit, direction)) {
				return { L, L, "No Linear Interval" };
			}
		}

		int final_M = curr_M;
		while (final_M > L) {
			if (_check_substrate_depletion(absorbances[final_M], depletion_limit, direction)) {
				final_M--;
			}
			else {
				break;
			}
		}
		return { L, final_M, "OK" };
	}

	std::string _linearity_check(const Vector& absorbances, int L, int M, double total_slope, double cycle_time, double linearity_limit = 0.2) {
		if (std::abs(total_slope) <= MIN_slope_CHECK_THRESHOLD) {
			return "Pass (Slope Low)";
		}
		int points_count = M - L + 1;
		std::vector<int> idx_f, idx_b;
		if (points_count > 8) {
			for (int i = 0; i < 6; ++i) idx_f.push_back(L + i);
			for (int i = 0; i < 6; ++i) idx_b.push_back(M - 5 + i);
		}
		else if (points_count >= 4) {
			for (int i = 0; i < 3; ++i) idx_f.push_back(L + i);
			for (int i = 0; i < 3; ++i) idx_b.push_back(M - 2 + i);
		}
		else {
			return "Skipped (N<4)";
		}

		double slope_f = _calculate_slope_least_squares(absorbances, idx_f, cycle_time);
		double slope_b = _calculate_slope_least_squares(absorbances, idx_b, cycle_time);

		double numerator = std::abs(slope_f - slope_b);
		if (numerator <= MIN_slope_CHECK_THRESHOLD) {
			return "Pass (Diff Low)";
		}

		double linearity = numerator / std::abs(total_slope);
		if (linearity > linearity_limit) {
			std::stringstream ss;
			ss << "Fail (LIN) val=" << std::fixed << std::setprecision(2) << linearity;
			return ss.str();
		}
		return "Pass";
	}

	std::pair<double, std::string> _enzyme_expansion(const Vector& absorbances, int delay_start_idx, int L, double cycle_time, double limit, int direction) {
		int search_end = L;
		while (search_end > delay_start_idx) {
			if (limit != -1.0 && _check_substrate_depletion(absorbances[search_end], limit, direction)) {
				search_end--;
			}
			else {
				break;
			}
		}
		if (search_end - delay_start_idx < 1) {
			return { 0.0, "ENC (Expansion Calc Fail)" };
		}

		double max_slope = 0.0;
		bool found = false;
		for (int i = delay_start_idx; i < search_end; ++i) {
			double slope = 60.0 * (absorbances[i + 1] - absorbances[i]) / cycle_time;
			if (!found || std::abs(slope) > std::abs(max_slope)) {
				max_slope = slope;
				found = true;
			}
		}
		return { max_slope, "EXP" };
	}


public:
	std::pair<double, std::string> calculate_response(const Vector& absorbances, int N, int P, int L, int M, double K, double cycle_time = 1.0,
		double depletion_limit = -1.0, int direction = 1, double lin_limit = 0.2) {

		LinearRangeResult range = _determine_linear_range(absorbances, L, M, depletion_limit, direction);
		if (range.status == "No Linear Interval") {
			return { 0.0, "Error: Substrate Depleted Early" };
		}

		int real_L = range.L;
		int real_M = range.M;
		int points_count = real_M - real_L + 1;

		if (points_count < 2) {
			int delay_start = (P > 0) ? P : 0;
			return _enzyme_expansion(absorbances, delay_start, L, cycle_time, depletion_limit, direction);
		}

		std::vector<int> reaction_indices;
		for (int i = real_L; i <= real_M; ++i) reaction_indices.push_back(i);
		double slope_reaction = _calculate_slope_least_squares(absorbances, reaction_indices, cycle_time);

		double slope_blank = 0.0;
		if (P > N) {
			std::vector<int> blank_indices;
			for (int i = N; i <= P; ++i) blank_indices.push_back(i);
			slope_blank = _calculate_slope_least_squares(absorbances, blank_indices, cycle_time);
		}

		double R = slope_reaction - (K * slope_blank);

		std::string lin_status = "Skipped";
		if (points_count >= 4) {
			lin_status = _linearity_check(absorbances, real_L, real_M, slope_reaction, cycle_time, lin_limit);
			if (lin_status.find("Fail") != std::string::npos) {
				std::stringstream ss;
				ss << "Result: " << std::fixed << std::setprecision(4) << R << " (Flag: LIN)";
				return { R, ss.str() };
			}
		}

		if (points_count == 2) {
			std::stringstream ss;
			ss << "Result: " << std::fixed << std::setprecision(4) << R << " (Flag: Points<3 Warning)";
			return { R, ss.str() };
		}

		std::stringstream ss;
		ss << "Result: " << std::fixed << std::setprecision(4) << R << " (Linearity: " << lin_status << ")";
		return { R, ss.str() };
	}
};

/**
 * @brief 存储速率法计算结果的结构体
 */
struct RateResult {
	double R;        ///< 最终反应度 (ΔAbs/min)
	double slope;    ///< 最小二乘法拟合的斜率 (Abs/sec)
	int n_points;    ///< 参与计算的点数
};


/**
 * @brief 速率法 (简化动力学法) 反应度计算模块
 * @details 通常指使用最小二乘法对指定区间的数据进行线性拟合来求速率。
 *          是 KineticCalculator 的一种简化形式。
 */
class RateMethodCalculator {
private:
	/**
 * @brief 计算速率法的反应度
 * @param[in] absorbances 吸光度数据列表
 * @param[in] start_idx 计算区间的起始点索引 (L)
 * @param[in] end_idx 计算区间的结束点索引 (M)
 * @param[in] cycle_time 测光周期(秒)
 * @return RateResult 包含反应度、斜率和点数的结构体
 */
	double _linear_regression(const Vector& x, const Vector& y) {
		size_t n = x.size();
		if (n < 2) return 0.0;
		double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
		double sum_y = std::accumulate(y.begin(), y.end(), 0.0);
		double sum_xy = 0.0;
		double sum_xx = 0.0;
		for (size_t i = 0; i < n; ++i) {
			sum_xy += x[i] * y[i];
			sum_xx += x[i] * x[i];
		}
		double denominator = n * sum_xx - sum_x * sum_x;
		return (denominator == 0) ? 0.0 : (n * sum_xy - sum_x * sum_y) / denominator;
	}

public:
	RateResult calculate_response(const Vector& absorbances, int start_idx, int end_idx, double cycle_time = 1.0) {
		if (start_idx >= absorbances.size() || end_idx >= absorbances.size() || start_idx > end_idx) {
			std::cerr << "[错误] 索引超出范围" << std::endl;
			return { 0.0, 0.0, 0 };
		}
		Vector y_data(absorbances.begin() + start_idx, absorbances.begin() + end_idx + 1);
		int n_points = y_data.size();

		if (n_points < 3) {
			std::cout << "[警告] 点数过少，建议使用至少3个点进行拟合" << std::endl;
			if (n_points == 2) {
				return { (y_data.back() - y_data.front()) / cycle_time * 60.0, 0.0, n_points };
			}
			return { 0.0, 0.0, 0 };
		}

		Vector x_data(n_points);
		for (int i = 0; i < n_points; ++i) x_data[i] = i * cycle_time;

		double slope = _linear_regression(x_data, y_data);
		return { slope * 60.0, slope, n_points };
	}


	/**
 * @brief (进阶功能) 检查反应曲线的线性度
 * @details 将计算区间分为前后两半，分别计算斜率并比较差异。
 * @param[in] absorbances 吸光度数据列表
 * @param[in] start_idx 检查区间的起始点索引
 * @param[in] end_idx 检查区间的结束点索引
 * @return bool 如果线性度在可接受范围内，返回 true，否则返回 false。
 */
	bool linearity_check(const Vector& absorbances, int start_idx, int end_idx) {
		int mid = (start_idx + end_idx) / 2;
		if (mid <= start_idx) return true;

		RateResult res1 = calculate_response(absorbances, start_idx, mid);
		RateResult res2 = calculate_response(absorbances, mid, end_idx);

		if (std::abs(res1.R) > 0.001) {
			double diff = std::abs(res1.R - res2.R) / std::abs(res1.R);
			if (diff > 0.1) {
				std::cout << "[非线性警告] 前段速率:" << std::fixed << std::setprecision(4) << res1.R
					<< ", 后段速率:" << res2.R << ", 偏差:" << std::setprecision(1) << diff * 100 << "%" << std::endl;
				return false;
			}
		}
		return true;
	}
};

/**
 * @brief 存储终点法计算结果的结构体
 */
struct EndpointResult {
	double R;  ///< 最终反应度
	double Ai; ///< 反应终点吸光度均值
	double Ab; ///< 空白读数吸光度均值
};


/**
 * @brief 终点法反应度计算模块
 * @details 对应文档 4.3.5 和 4.3.6。
 *          适用于反应达到或接近终点时进行测量的项目。
 */
class EndpointCalculator {
private:
	/**
 * @brief 计算基础反应度 R
 * @details 公式为: R = Ai - k * Ab
 * @param[in] absorbances 反应曲线数据
 * @param[in] N 空白时间起始点索引
 * @param[in] P 空白时间结束点索引
 * @param[in] L 反应时间起始点索引
 * @param[in] M 反应时间结束点索引
 * @param[in] K 体积校正因子 (k)
 * @return EndpointResult 包含最终反应度 R 和中间值的结构体
 */
	double _calculate_mean_abs(const Vector& absorbances, int start_idx, int end_idx) {
		if (absorbances.empty() || start_idx >= absorbances.size()) return 0.0;
		int real_end = std::min((size_t)end_idx, absorbances.size() - 1);
		if (start_idx > real_end) return 0.0;

		Vector segment(absorbances.begin() + start_idx, absorbances.begin() + real_end + 1);
		return std::accumulate(segment.begin(), segment.end(), 0.0) / segment.size();
	}
public:
	EndpointResult calculate_raw_response(const Vector& absorbances, int N, int P, int L, int M, double K) {
		double Ai = _calculate_mean_abs(absorbances, L, M);
		double Ab = (N == 0 && P == 0 && K == 0) ? 0.0 : _calculate_mean_abs(absorbances, N, P);
		double R = Ai - (K * Ab);
		return { R, Ai, Ab };
	}
	/**
 * @brief 计算样本空白校正后的反应度 R'
 * @details 公式为: R' = R_sample - R_sample_blank
 * @param[in] sample_data 样本反应曲线
 * @param[in] sample_blank_data 样本空白反应曲线 (如果无，则为空 vector)
 * @param[in] N, P, L, M, K 同 `calculate_raw_response` 的计算参数
 * @return double 校正后的最终反应度 R'
 */
	double calculate_corrected_response(const Vector& sample_data, const Vector& sample_blank_data, int N, int P, int L, int M, double K) {
		EndpointResult res_sample = calculate_raw_response(sample_data, N, P, L, M, K);
		if (sample_blank_data.empty()) return res_sample.R;

		EndpointResult res_sb = calculate_raw_response(sample_blank_data, N, P, L, M, K);
		double R_prime = res_sample.R - res_sb.R;

		std::cout << "[计算细节] R样本=" << std::fixed << std::setprecision(4) << res_sample.R
			<< ", R样空=" << res_sb.R << ", R'=" << R_prime << std::endl;
		return R_prime;
	}
};


/**
 * @brief 线性校准模块
 * @details 对应文档 4.6.1。
 *          通过标准品建立反应度(R)与浓度(C)之间的线性关系 C = K * (R - R0)。
 *          支持单点、两点和多点线性校准。
 */
class LinearCalibrator {
private:
	double K = 1.0;
	double R0 = 0.0;
	bool is_calibrated = false;
	const double SCALE_FACTOR = 10000.0;

	double _scale_R(double raw_R) { return raw_R / SCALE_FACTOR; }

public:
	/**
 * @brief 单点线性校准 (K因数法)
 * @param[in] user_K 用户输入的 K 因数
 * @param[in] user_R0_raw 试剂空白的原始反应度 (未除以10000)
 */
	void fit_single_point(double user_K, double user_R0_raw = 0.0) {
		this->K = user_K;
		this->R0 = _scale_R(user_R0_raw);
		this->is_calibrated = true;
		std::cout << "[校准成功] 单点法: K=" << std::fixed << std::setprecision(4) << this->K << ", R0=" << this->R0 << std::endl;
	}
	/**
 * @brief 两点线性校准
 * @param[in] concs 包含两个校准品浓度的向量
 * @param[in] resps_raw 包含两个校准品原始反应度的向量
 */
	void fit_two_point(const Vector& concs, const Vector& resps_raw) {
		if (concs.size() != 2 || resps_raw.size() != 2) {
			std::cerr << "[错误] 两点法必须传入 2 个点的浓度和反应度" << std::endl;
			return;
		}
		double R1 = _scale_R(resps_raw[0]), R2 = _scale_R(resps_raw[1]);
		if (std::abs(R2 - R1) < 1e-9) {
			std::cerr << "[错误] 两个校准品的反应度相同，无法计算斜率" << std::endl;
			return;
		}
		this->K = (concs[1] - concs[0]) / (R2 - R1);
		this->R0 = (this->K != 0) ? (R1 - (concs[0] / this->K)) : 0.0;
		this->is_calibrated = true;
		std::cout << "[校准成功] 两点法: K=" << std::fixed << std::setprecision(4) << this->K << ", R0=" << this->R0 << std::endl;
	}
	/**
 * @brief 多点线性校准 (最小二乘法)
 * @param[in] concs 包含多个校准品浓度的向量 (n >= 3)
 * @param[in] resps_raw 包含多个校准品原始反应度的向量
 */
	void fit_multi_point(const Vector& concs, const Vector& resps_raw) {
		size_t n = concs.size();
		if (n < 3 || concs.size() != resps_raw.size()) {
			std::cerr << "[错误] 多点法至少需要 3 个校准品" << std::endl;
			return;
		}
		Vector R;
		for (double r : resps_raw) R.push_back(_scale_R(r));

		double sum_C = std::accumulate(concs.begin(), concs.end(), 0.0);
		double sum_R = std::accumulate(R.begin(), R.end(), 0.0);
		double sum_CR = 0.0, sum_R2 = 0.0;
		for (size_t i = 0; i < n; ++i) {
			sum_CR += concs[i] * R[i];
			sum_R2 += R[i] * R[i];
		}

		double denominator = sum_R2 - (sum_R * sum_R) / n;
		if (std::abs(denominator) < 1e-9) {
			std::cerr << "[错误] 反应度数据方差为0，无法拟合" << std::endl;
			return;
		}

		double numerator = sum_CR - (sum_C * sum_R) / n;
		this->K = numerator / denominator;
		this->R0 = (this->K != 0) ? (sum_R / n) - (sum_C / n) / this->K : 0.0;
		this->is_calibrated = true;
		std::cout << "[校准成功] 多点法(" << n << "点): K=" << std::fixed << std::setprecision(4) << this->K << ", R0=" << this->R0 << std::endl;
	}
	/**
 * @brief 使用已建立的线性模型计算样本浓度
 * @note 输入的样本反应度 `sample_resp_raw` 是未缩放的原始值。
 * @param[in] sample_resp_raw 样本的原始反应度
 * @return double 计算出的样本浓度
 */
	double calculate_concentration(double sample_resp_raw) {
		if (!is_calibrated) {
			std::cerr << "[警告] 尚未校准" << std::endl;
			return 0.0;
		}
		double R = _scale_R(sample_resp_raw);
		return this->K * (R - this->R0);
	}
};


int main() {
	// 设置输出格式
	std::cout << std::fixed << std::setprecision(4);

	// ========================================================
	// 模块一: 固定时间法 (Fixed-Time Method)
	// ========================================================
	std::cout << "\n" << std::string(50, '=') << std::endl;
	std::cout << "模块一: 固定时间法 (Fixed-Time Method)" << std::endl;
	std::cout << std::string(50, '=') << std::endl;
	std::cout << "适用于反应速率在监测时间内恒定的项目，如肌酐(苦味酸法)。\n" << std::endl;

	double CYCLE_TIME_SEC = 18.0;
	FixedTimeCalculator calc_ft;
	Vector raw_curve_ft;
	for (int i = 0; i < 17; ++i) raw_curve_ft.push_back(0.10 + i * 0.001);
	for (int i = 0; i < 20; ++i) raw_curve_ft.push_back(0.12 + i * 0.02);

	std::cout << "-> 模拟数据: 肌酐反应曲线 (Jaffe Method)" << std::endl;
	std::cout << "  - 数据点总数: " << raw_curve_ft.size() << std::endl;
	std::cout << "  - 测光周期: " << CYCLE_TIME_SEC << " 秒\n" << std::endl;

	std::cout << "--- 场景 1: 双试剂，双区间扣除 (带体积修正) ---" << std::endl;
	std::cout << "   目的: 扣除主反应前由试剂1引起的非特异性反应速率。" << std::endl;
	int N1 = 5, P1 = 16, L1 = 17, M1 = 30;
	double K1 = 0.8;
	FixedTimeResult res1 = calc_ft.calculate_response(raw_curve_ft, N1, P1, L1, M1, K1, CYCLE_TIME_SEC);
	std::cout << "   输入参数:" << std::endl;
	std::cout << "     - 空白监测窗 (N-P): " << N1 << "-" << P1 << std::endl;
	std::cout << "     - 反应监测窗 (L-M): " << L1 << "-" << M1 << std::endl;
	std::cout << "     - 体积修正因子 (K): " << K1 << std::endl;
	std::cout << "   计算过程:" << std::endl;
	std::cout << "     - 反应区间斜率: " << std::setprecision(6) << res1.slope_reaction << " Abs/sec" << std::endl;
	std::cout << "     - 空白区间斜率: " << res1.slope_blank << " Abs/sec" << std::endl;
	std::cout << "     - 计算公式: R = 60 * (反应斜率 - K * 空白斜率)" << std::endl;
	std::cout << "     -           = 60 * (" << res1.slope_reaction << " - " << K1 << " * " << res1.slope_blank << ")" << std::endl;
	std::cout << "   最终结果:" << std::endl;
	std::cout << "     - 最终反应度 R: " << std::setprecision(4) << res1.R << " ΔAbs/min\n" << std::endl;

	std::cout << "--- 场景 2: 不扣除空白 (K=0) ---" << std::endl;
	std::cout << "   目的: 仅计算主反应区间的速率，适用于无明显试剂空白的项目。" << std::endl;
	FixedTimeResult res2 = calc_ft.calculate_response(raw_curve_ft, 0, 0, 17, 30, 0.0, CYCLE_TIME_SEC);
	std::cout << "   输入参数:" << std::endl;
	std::cout << "     - 空白监测窗 (N-P): 0-0 (不使用)" << std::endl;
	std::cout << "     - 反应监测窗 (L-M): 17-30" << std::endl;
	std::cout << "     - 体积修正因子 (K): 0.0" << std::endl;
	std::cout << "   计算过程:" << std::endl;
	std::cout << "     - 反应区间斜率: " << std::setprecision(6) << res2.slope_reaction << " Abs/sec" << std::endl;
	std::cout << "     - 计算公式: R = 60 * 反应斜率" << std::endl;
	std::cout << "   最终结果:" << std::endl;
	std::cout << "     - 最终反应度 R: " << std::setprecision(4) << res2.R << " ΔAbs/min\n" << std::endl;

	// ========================================================
	// 模块二: 终点法 (Endpoint Method)
	// ========================================================
	std::cout << "\n" << std::string(50, '=') << std::endl;
	std::cout << "模块二: 终点法 (Endpoint Method)" << std::endl;
	std::cout << std::string(50, '=') << std::endl;
	std::cout << "适用于反应达到或接近终点时进行测量的项目，如总蛋白、白蛋白等。\n" << std::endl;

	EndpointCalculator calc_ep;
	Vector raw_curve_ep;
	for (int i = 0; i < 10; ++i) raw_curve_ep.push_back(0.1);
	for (int i = 0; i < 30; ++i) raw_curve_ep.push_back(0.12);
	for (int i = 0; i < 30; ++i) raw_curve_ep.push_back(0.5 + 0.01 * i);
	std::cout << "-> 模拟数据: 一个典型的终点法反应曲线" << std::endl;
	std::cout << "  - 数据点总数: " << raw_curve_ep.size() << "\n" << std::endl;

	std::cout << "--- 场景 1: 双试剂，空白时间在反应启动前 (带体积修正) ---" << std::endl;
	EndpointResult ep_res1 = calc_ep.calculate_raw_response(raw_curve_ep, 5, 15, 50, 60, 0.8);
	std::cout << "   输入参数:" << std::endl;
	std::cout << "     - 空白窗 (N-P): 5-15" << std::endl;
	std::cout << "     - 反应窗 (L-M): 50-60" << std::endl;
	std::cout << "     - 体积修正因子 (K): 0.8" << std::endl;
	std::cout << "   计算过程:" << std::endl;
	std::cout << "     - 反应终点均值 (Ai): " << ep_res1.Ai << std::endl;
	std::cout << "     - 空白读数均值 (Ab): " << ep_res1.Ab << std::endl;
	std::cout << "     - 计算公式: R = Ai - K * Ab" << std::endl;
	std::cout << "   最终结果:" << std::endl;
	std::cout << "     - 最终反应度 R = " << ep_res1.Ai << " - 0.8 * " << ep_res1.Ab << " = " << ep_res1.R << "\n" << std::endl;

	std::cout << "--- 场景 4: 样本空白校正 ---" << std::endl;
	std::cout << "   目的: 扣除由样本自身颜色、浊度等引起的背景吸收。" << std::endl;
	Vector sample_curve;
	for (double x : raw_curve_ep) sample_curve.push_back(x + 0.5);
	Vector sample_blank_curve(70, 0.52);
	double final_R = calc_ep.calculate_corrected_response(sample_curve, sample_blank_curve, 5, 15, 50, 60, 0.8);
	std::cout << "   说明: 使用场景1的参数，对一个背景增高(如脂血)的样本进行校正。" << std::endl;
	std::cout << "   最终结果:" << std::endl;
	std::cout << "     - 校正后反应度 R': " << final_R << "\n" << std::endl;

	// ========================================================
	// 模块三: 动力学法 (Kinetic Method)
	// ========================================================
	std::cout << "\n" << std::string(50, '=') << std::endl;
	std::cout << "模块三: 动力学法 (Kinetic Method)" << std::endl;
	std::cout << std::string(50, '=') << std::endl;
	std::cout << "主要用于酶活性检测，通过监测反应速率来计算。包含复杂的线性检查和范围调整逻辑。\n" << std::endl;

	KineticCalculator calc_kin;
	double CYCLE_TIME_KIN = 10.0;
	Vector raw_data_kin = { 2.5, 2.2, 1.8, 1.4, 1.0, 0.6, 0.4, 0.3, 0.28, 0.28, 0.28, 0.28, 0.28 };
	double LIMIT = 0.5;
	int DIR = -1;
	std::cout << "-> 模拟数据: 高活性酶反应，导致底物耗尽" << std::endl;
	std::cout << "  - 反应方向: 下降 (DIR=" << DIR << ")" << std::endl;
	std::cout << "  - 耗尽限值: " << LIMIT << " Abs" << std::endl;
	std::cout << "  - 测光周期: " << CYCLE_TIME_KIN << " 秒\n" << std::endl;

	std::cout << "--- 测试 1: 底物耗尽，触发线性范围缩减 ---" << std::endl;
	auto res_kin1 = calc_kin.calculate_response(raw_data_kin, 0, 0, 3, 10, 0, CYCLE_TIME_KIN, LIMIT, DIR);
	std::cout << "   输入参数:" << std::endl;
	std::cout << "     - 预设区间 (L-M): [3~10]" << std::endl;
	std::cout << "   运行信息: " << res_kin1.second << std::endl;
	std::cout << "   说明: 算法检测到在 M=10 点时反应已耗尽，自动缩短计算区间以获得有效速率。\n" << std::endl;

	std::cout << "--- 测试 2: 起点耗尽，触发酶线性扩展 (EXP) ---" << std::endl;
	auto res_kin2 = calc_kin.calculate_response(raw_data_kin, 0, 0, 6, 12, 0, CYCLE_TIME_KIN, LIMIT, DIR);
	std::cout << "   输入参数:" << std::endl;
	std::cout << "     - 预设区间 (L-M): [6~12]" << std::endl;
	std::cout << "   运行信息: " << res_kin2.second << std::endl;
	std::cout << "   说明: 算法检测到在 L=6 点时反应就已耗尽，预设主区间完全无效。" << std::endl;
	std::cout << "           自动回溯到延迟期，寻找最快的反应速率作为结果，并标记为 'EXP'。\n" << std::endl;

	// ========================================================
	// 模块四: 速率法 (Rate Method)
	// ========================================================
	std::cout << "\n" << std::string(50, '=') << std::endl;
	std::cout << "模块四: 速率法 (Rate Method)" << std::endl;
	std::cout << std::string(50, '=') << std::endl;
	std::cout << "是动力学法的一种简化形式，通常指用最小二乘法对指定区间的数据进行线性拟合来求速率。\n" << std::endl;

	RateMethodCalculator calc_rate;
	Vector raw_curve_rate;
	for (int i = 0; i < 50; ++i) raw_curve_rate.push_back(1.5 + (-0.002 * i * 18.0) + (rand() % 1000 / 1000000.0 - 0.0005));
	std::cout << "-> 模拟数据: ALT (谷丙转氨酶) 反应曲线 (带随机噪声)" << std::endl;
	std::cout << "  - 测光周期: 18.0 秒" << std::endl;
	std::cout << "  - 理论斜率: -0.002 Abs/sec\n" << std::endl;

	std::cout << "--- 速率法计算演示 ---" << std::endl;
	RateResult rate_res = calc_rate.calculate_response(raw_curve_rate, 10, 30, 18.0);
	bool is_linear = calc_rate.linearity_check(raw_curve_rate, 10, 30);
	std::cout << "   输入参数:" << std::endl;
	std::cout << "     - 计算区间 (L-M): 10-30" << std::endl;
	std::cout << "   计算结果:" << std::endl;
	std::cout << "     - 参与拟合点数: " << rate_res.n_points << std::endl;
	std::cout << "     - 最小二乘法拟合斜率: " << std::setprecision(6) << rate_res.slope << " Abs/sec" << std::endl;
	std::cout << "     - 最终反应度 R: " << std::setprecision(4) << rate_res.R << " ΔAbs/min" << std::endl;
	std::cout << "   对比:" << std::endl;
	std::cout << "     - 理论真值: " << (-0.002 * 60.0) << " ΔAbs/min" << std::endl;
	std::cout << "   检查项:" << std::endl;
	std::cout << "     - 线性检查结果: " << (is_linear ? "通过" : "失败") << "\n" << std::endl;

	// ========================================================
	// 模块五: 线性校准 (Linear Calibration)
	// ========================================================
	std::cout << "\n" << std::string(50, '=') << std::endl;
	std::cout << "模块五: 线性校准 (Linear Calibration)" << std::endl;
	std::cout << std::string(50, '=') << std::endl;
	std::cout << "通过标准品建立反应度(R)与浓度(C)之间的线性关系 C = a*R + b。\n" << std::endl;

	LinearCalibrator calibrator;
	std::cout << "--- 测试 1: 单点线性校准 (K因数法) ---" << std::endl;
	calibrator.fit_single_point(100.0, 500);
	std::cout << "   校准参数: K = 100.0, R0 = 500 (原始信号值)" << std::endl;
	std::cout << "   样本计算: 原始反应度 2000 -> 计算浓度: " << calibrator.calculate_concentration(2000) << "\n" << std::endl;

	std::cout << "--- 测试 2: 两点线性校准 ---" << std::endl;
	calibrator.fit_two_point({ 10.0, 50.0 }, { 1000, 5000 });
	std::cout << "   校准点: C1=10.0 @ R1=1000, C2=50.0 @ R2=5000" << std::endl;
	std::cout << "   样本计算: 原始反应度 2000 -> 计算浓度: " << calibrator.calculate_concentration(2000) << "\n" << std::endl;

	std::cout << "--- 测试 3: 多点线性校准 (最小二乘法) ---" << std::endl;
	calibrator.fit_multi_point({ 10.0, 30.0, 50.0, 68.0 }, { 1000, 3000, 5000, 7000 });
	std::cout << "   校准点 (带噪声): C={10.0, 30.0, 50.0, 68.0} @ R={1000, 3000, 5000, 7000}" << std::endl;
	std::cout << "   样本计算: 原始反应度 4000 -> 计算浓度: " << calibrator.calculate_concentration(4000) << "\n" << std::endl;

	// ========================================================
	// 模块六: 非线性校准 (Non-Linear Calibration)
	// ========================================================
	std::cout << "\n" << std::string(50, '=') << std::endl;
	std::cout << "模块六: 非线性校准 (Non-Linear Calibration)" << std::endl;
	std::cout << std::string(50, '=') << std::endl;
	std::cout << "用于免疫等反应，其反应度(R)与浓度(C)不成线性关系，需要用复杂曲线拟合。\n" << std::endl;
	NonLinearCalibrator calib_nl;

	std::cout << "--- 测试 1: Polynomial 5P 模型 ---" << std::endl;
	calib_nl.fit_polynomial_5p({ 0, 10, 50, 100, 200 }, { 500, 1200, 2500, 3200, 3800 });
	std::cout << "   校准点: C={0, 10, 50, 100, 200} @ R={500, 1200, 2500, 3200, 3800}" << std::endl;
	std::cout << "   样本计算: R=2000 -> Conc=" << calib_nl.calculate_concentration(2000) << " (预期在 10-50 之间)\n" << std::endl;

	std::cout << "--- 测试 2: Parabola (二次多项式) 模型 ---" << std::endl;
	calib_nl.fit_parabola({ 0, 50, 100 }, { 100, 5000, 9000 });
	std::cout << "   校准点: C={0, 50, 100} @ R={100, 5000, 9000}" << std::endl;
	std::cout << "   样本计算: R=3000 -> Conc=" << calib_nl.calculate_concentration(3000) << "\n" << std::endl;

	std::cout << "--- 测试 3: Logistic 4P 模型 (直接使用参数) ---" << std::endl;
	calib_nl.set_params_logistic_4p(100, 5000, 2.0, -1.0);
	std::cout << "   模型参数: R0=100, K=5000, a=2.0, b=-1.0" << std::endl;
	std::cout << "   样本计算: R=2600 -> Conc=" << calib_nl.calculate_concentration(2600) << "\n" << std::endl;

	return 0;
}
