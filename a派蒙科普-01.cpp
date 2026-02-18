//算法由AI提供
//仅供测试，严禁用于生产环境
//本文件在CC-BY-4.0许可协议下分发:https://creativecommons.org/licenses/by/4.0/legalcode.zh-hans
#include <iostream>
#include <cmath>
#include <unordered_map>
#include <limits>
#include <stdexcept>
#include <ctime>

// 高精度数，以 2 为底的科学计数法表示： mantissa * 2^exponent
struct HighPrecision {
    bool sign;          // true 为正，false 为负（本实现中只用正数）
    double mantissa;    // 尾数，范围 [1, 2)
    long long exponent; // 指数，可为任意整数（包括很大的负数）

    HighPrecision() : sign(true), mantissa(1.0), exponent(-std::numeric_limits<long long>::max()) {}

    HighPrecision(double val) {
        if (val == 0.0) {
            mantissa = 1.0;
            exponent = -std::numeric_limits<long long>::max();
            sign = true;
            return;
        }
        sign = val > 0;
        val = std::fabs(val);
        int exp;
        mantissa = std::frexp(val, &exp);   // 返回 [0.5, 1)
        exponent = exp - 1;                  // 调整到 [1, 2)
        mantissa *= 2;
    }

    HighPrecision(double m, long long e, bool s = true) : sign(s), mantissa(m), exponent(e) {
        normalize();
    }

    void normalize() {
        if (mantissa == 0.0 || exponent == -std::numeric_limits<long long>::max()) {
            exponent = -std::numeric_limits<long long>::max();
            mantissa = 1.0;
            return;
        }
        while (mantissa >= 2.0) {
            mantissa *= 0.5;
            exponent++;
        }
        while (mantissa < 1.0 && exponent > -std::numeric_limits<long long>::max()) {
            mantissa *= 2.0;
            exponent--;
        }
    }

    HighPrecision mul_half() const {
        HighPrecision res = *this;
        res.exponent--;
        return res;
    }

    HighPrecision add(const HighPrecision& other) const {
        if (exponent == -std::numeric_limits<long long>::max()) return other;
        if (other.exponent == -std::numeric_limits<long long>::max()) return *this;

        long long target_exp = std::max(exponent, other.exponent);
        double m1 = mantissa * std::ldexp(1.0, exponent - target_exp);
        double m2 = other.mantissa * std::ldexp(1.0, other.exponent - target_exp);
        double sum = m1 + m2;
        return HighPrecision(sum, target_exp);
    }

    HighPrecision sub(const HighPrecision& other) const {
        if (other.exponent == -std::numeric_limits<long long>::max()) return *this;
        if (exponent == -std::numeric_limits<long long>::max())
            throw std::runtime_error("Negative result in sub");

        long long target_exp = std::max(exponent, other.exponent);
        double m1 = mantissa * std::ldexp(1.0, exponent - target_exp);
        double m2 = other.mantissa * std::ldexp(1.0, other.exponent - target_exp);
        double diff = m1 - m2;
        if (diff < 0) throw std::runtime_error("Negative diff in sub");
        return HighPrecision(diff, target_exp);
    }

    int cmp(const HighPrecision& other) const {
        if (exponent > other.exponent) return 1;
        if (exponent < other.exponent) return -1;
        if (mantissa > other.mantissa) return 1;
        if (mantissa < other.mantissa) return -1;
        return 0;
    }

    bool is_zero() const {
        return exponent == -std::numeric_limits<long long>::max();
    }

    double to_double() const {
        if (exponent == -std::numeric_limits<long long>::max()) return 0.0;
        return std::ldexp(mantissa, exponent);
    }

    void print() const {
        if (exponent == -std::numeric_limits<long long>::max()) {
            std::cout << "0";
            return;
        }
        std::cout << mantissa << " * 2^" << exponent;
    }
};

std::unordered_map<double, HighPrecision> cache;

int leadingOnes(double y) {
    int k = 0;
    while (y >= 0.5) {
        y = 2 * y - 1;
        k++;
    }
    return k;
}

// f(x) for x in [0,1)
HighPrecision f0(double x) {
    int k = leadingOnes(x);
    double val = 1.0 - x - std::ldexp(1.0, -k - 1);
    return HighPrecision(val);
}

// f(x) for x in [1,2)   x = 1 + y, y in [0,1)
HighPrecision f1(double y) {
    int m = 0;
    double cur = y;
    while (true) {
        int nu = leadingOnes(cur);
        double S = 2 * cur + std::ldexp(1.0, -nu - 1);
        if (S < 1.0) {
            HighPrecision res = f0(S);
            for (int i = 0; i < m + 1; ++i)
                res = res.mul_half();
            return res;
        }
        else {
            cur = S - 1.0;
            m++;
        }
    }
}

// f(x) for x in [2,3)   x = 2 + y, y in [0,1)
HighPrecision f2(double y, int max_steps = 1000000) {
    int m = 0;
    HighPrecision cur_y = HighPrecision(y);
    while (m <= max_steps) {
        double cur_d = cur_y.to_double();
        HighPrecision h = f1(cur_d); // f(1+cur_y)
        if (h.is_zero()) {
            // h is extremely small, cannot reduce further
            return HighPrecision(0.0);
        }
        if (cur_y.cmp(h) >= 0) { // cur_y >= h
            HighPrecision new_cur = cur_y.sub(h);
            cur_y = new_cur;
            m++;
        }
        else {
            // cur_y < h, exit to [1,2)
            // x = 2 + cur_y, after subtraction: x' = 2 + cur_y - h = 1 + (1 + cur_y - h)
            // let u = 1 + cur_y - h ∈ (0,1). Then x' = 1 + u, need f1(u)
            double u = 1.0 + (cur_y.to_double() - h.to_double());
            HighPrecision res = f1(u);
            for (int i = 0; i < m + 1; ++i)
                res = res.mul_half();
            return res;
        }
    }
    return HighPrecision(0.0);
}

// f(x) for x in [3,4)   x = 3 + y, y in [0,1)
HighPrecision f3(double y, int max_steps = 1000000) {
    int m = 0;
    HighPrecision cur_y = HighPrecision(y);
    while (m <= max_steps) {
        double cur_d = cur_y.to_double();
        HighPrecision h = f2(cur_d); // f(2+cur_y)
        if (h.is_zero()) {
            return HighPrecision(0.0);
        }
        if (cur_y.cmp(h) >= 0) {
            HighPrecision new_cur = cur_y.sub(h);
            cur_y = new_cur;
            m++;
        }
        else {
            // cur_y < h, exit to [2,3)
            // x = 3 + cur_y, after subtraction: x' = 3 + cur_y - h = 2 + (1 + cur_y - h)
            // let v = 1 + cur_y - h ∈ (0,1). Then x' = 2 + v, need f2(v)
            double v = 1.0 + (cur_y.to_double() - h.to_double());
            HighPrecision res = f2(v);
            for (int i = 0; i < m + 1; ++i)
                res = res.mul_half();
            return res;
        }
    }
    return HighPrecision(0.0);
}

// f(x) for x in [4,4.5]   x = 4 + y, y in [0,0.5]
HighPrecision f4(double y, int max_steps = 1000000) {
    int m = 0;
    HighPrecision cur_y = HighPrecision(y);
    while (m <= max_steps) {
        double cur_d = cur_y.to_double();
        HighPrecision h = f3(cur_d); // f(3+cur_y)
        if (h.is_zero()) {
            return HighPrecision(0.0);
        }
        if (cur_y.cmp(h) >= 0) {
            HighPrecision new_cur = cur_y.sub(h);
            cur_y = new_cur;
            m++;
        }
        else {
            // cur_y < h, exit to [3,4)
            // x = 4 + cur_y, after subtraction: x' = 4 + cur_y - h = 3 + (1 + cur_y - h)
            // let w = 1 + cur_y - h ∈ (0,1). Then x' = 3 + w, need f3(w)
            double w = 1.0 + (cur_y.to_double() - h.to_double());
            HighPrecision res = f3(w);
            for (int i = 0; i < m + 1; ++i)
                res = res.mul_half();
            return res;
        }
    }
    return HighPrecision(0.0);
}

HighPrecision f(double x) {
    if (x < 0) {
        return HighPrecision(-x);
    }
    auto it = cache.find(x);
    if (it != cache.end()) return it->second;

    HighPrecision result;
    if (x < 1.0) {
        result = f0(x);
    }
    else if (x < 2.0) {
        result = f1(x - 1.0);
    }
    else if (x < 3.0) {
        result = f2(x - 2.0);
    }
    /*else if (x < 4.0) {
        result = f3(x - 3.0);
    }
    else if (x <= 4.5) {
        result = f4(x - 4.0);
    }*///超出3.0时，100%引发异常:被减数小于减数(cmp 误判),原因是值过于微小
    else {
        std::cout << "超出可接受范围\n";
        result = HighPrecision(0.0);
    }
    cache[x] = result;
    return result;
}

#include <chrono>

int main() {
    double x;
    std::cout << "Enter x (<3.0) repeatedly, Ctrl+C to exit." << std::endl;
    while (true) {
        std::cout << "x: ";
        if (!(std::cin >> x)) {
            break;
        }

        auto start = std::chrono::high_resolution_clock::now();
        HighPrecision val = f(x);
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::time_t t = std::time(nullptr);
        char time_str[100];
        std::strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S", std::localtime(&t));

        std::cout << "[" << time_str << "] f(" << x << ") = ";
        val.print();
        std::cout << " (耗时 " << duration.count() << " μs)" << std::endl;
        std::cout << "当最大值高于此值时，可能被截断：1 * 2^" << LLONG_MAX << std::endl;
        std::cout << "当最小值低于此值时，可能被截断：1 * 2^" << LLONG_MIN << std::endl;
    }
    return 0;
}
