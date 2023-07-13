#include "polynomials_field.h"

int polynomials_field::get(int i, const std::vector<int> &v) {
    return (0 <= i && i < v.size()) ? v[i] : 0;
}

void polynomials_field::shrink(std::vector<int> &a) {
    int erasable = 0;
    while (!a.empty() && erasable < a.size() && a[a.size() - 1 - erasable] == 0) {
        erasable++;
    }
    a.resize(a.size() - erasable);
}

[[nodiscard]] std::vector<int> polynomials_field::mul(const std::vector<int> &P, const std::vector<int> &Q) const {
    std::vector<int> ans(P.size() + Q.size() - 1);
    for (int i = 0; i < ans.size(); ++i) {
        int temp = 0;
        for (int j = 0; j <= i; ++j) {
            int element_p = get(j, P);
            int element_q = get(i - j, Q);
            if (element_p != 0 && element_q != 0) {
                int log_p = GF.getLog(element_p);
                int log_q = GF.getLog(element_q);
                if (log_p + log_q >= GF.q - 1) {
                    temp ^= GF.log_to_element[log_p + log_q - GF.q + 1];
                } else {
                    temp ^= GF.log_to_element[log_q + log_p];
                }
            }
        }
        ans[i] = temp;
    }
    return ans;
}

[[nodiscard]] std::vector<int>
polynomials_field::sum(const std::vector<int> &P, const std::vector<int> &Q) const {
    std::vector<int> ans(std::max(P.size(), Q.size()));
    for (int i = 0; i < ans.size(); ++i) {
        ans[i] = get(i, P) ^ get(i, Q);
    }
    return ans;
}

void polynomials_field::scalar_mul_vector(int scalar, std::vector<int> &vector, std::vector<int> &result) const {
    result.resize(vector.size());
    if (scalar == 0) {
        result.clear();
        return;
    }
    int log = GF.getLog(scalar);
    for (int i = 0; i < result.size(); ++i) {
        if (vector[i] == 0) {
            result[i] = 0;
            continue;
        }
        int elem_log = GF.getLog(vector[i]);
        if (elem_log + log < GF.q - 1) {
            result[i] = GF.log_to_element[log + elem_log];
        } else {
            result[i] = GF.log_to_element[log + elem_log - GF.q + 1];
        }
    }
}


std::vector<int> polynomials_field::shiftLeft(const std::vector<int> &a, size_t shift) {
    std::vector<int> result(a.size() + shift, 0);
    for (int i = 0; i < a.size(); ++i) {
        result[i + shift] = a[i];
    }
    return result;
}

std::vector<int> polynomials_field::mod(std::vector<int> a, std::vector<int> &b) const {
    std::vector<int> t;
    while (a.size() >= b.size()) {
        t = shiftLeft(b, a.size() - b.size());
        int coeff = GF.multiply_elements(a[a.size() - 1], GF.getInverse(t[t.size() - 1]));
        scalar_mul_vector(coeff, t, t);
        t = sum(a, t);
        shrink(t);
        a = t;
    }

    return a;
}

std::vector<int> polynomials_field::div(std::vector<int> a, std::vector<int> &b) const {
    shrink(a);
    shrink(b);
    if (a.size() < b.size()) {
        return {0};
    }
    std::vector<int> result(a.size() - b.size() + 1, 0);
    std::vector<int> t;
    while (!a.empty()) {
        t = shiftLeft(b, a.size() - b.size());
        int coeff = GF.multiply_elements(a[a.size() - 1], GF.getInverse(t[t.size() - 1]));
        result[a.size() - b.size()] = coeff;

        t = mul({coeff}, t);
        a = sum(t, a);
        shrink(a);
    }
    return result;
}

int polynomials_field::degree(std::vector<int> &polynomial) {
    if (polynomial.empty()) {
        return -INF;
    }
    for (int i = polynomial.size() - 1; i >= 0; --i) {
        if (polynomial[i] != 0) {
            return i;
        }
    }
    return -INF;

}

std::vector<int> polynomials_field::gcd(std::vector<int> a, std::vector<int> b) const {
    shrink(a);
    shrink(b);

    while (!b.empty()) {
        a = mod(a, b);
        std::swap(a, b);
    }
    return a;
}

int polynomials_field::evaluate(const std::vector<int> &polynomial, int a) const {
    int value = 0;
    if (a == 0) {
        return polynomial[0];
    }

    int power = GF.element_to_log[a];
    int degree = 0;
    for (int j = 0, log_j = 0; j < polynomial.size(); ++j) {
        if (polynomial[j] != 0) {
            int log = GF.element_to_log[polynomial[j]];
            if (log + degree < GF.q - 1) {
                value ^= GF.log_to_element[log + degree];
            } else {
                value ^= GF.log_to_element[log + degree - GF.q + 1];
            }
        }
        log_j++;
        degree += power;
        if (degree >= GF.q - 1) {
            degree -= GF.q - 1;
        }
        if (log_j == GF.q - 1) {
            log_j = 0;
            degree = 0;
        }

    }


    return value;
}

std::vector<int> polynomials_field::get_roots(std::vector<int> &polynomial) const {
    std::vector<int> roots(polynomial.size() - 1, GF.q);
    int roots_count = 0;
    size_t polynomial_degree = polynomial.size() - 1;
    for (int i = 0; i < GF.q; ++i) {
        if (evaluate(polynomial, i) == 0) {
            roots[roots_count] = i;
            roots_count++;
            if (roots_count == polynomial_degree) {
                break;
            }
        }
    }
    roots.resize(roots_count);
    return roots;
}
