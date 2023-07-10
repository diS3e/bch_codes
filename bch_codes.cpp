//
// Created by s3 on 23.04.23.
//
#include "bch_codes.h"
#include <algorithm>
#include <set>
#include <stdexcept>
#include <iostream>

int bch_code::get(int i, const std::vector<int> &v) {
    return (0 <= i && i < v.size()) ? v[i] : 0;
}

// Метод shrink - убираются нулевые старшие степени
void bch_code::shrink(std::vector<int> &a) {
    int erasable = 0;
    while (!a.empty() && erasable < a.size() && a[a.size() - 1 - erasable] == 0) {
        erasable++;
    }
    a.resize(a.size() - erasable);
}

[[nodiscard]] std::vector<int>
bch_code::multiply_polynomial(const std::vector<int> &P, const std::vector<int> &Q) const {
    std::vector<int> ans(P.size() + Q.size() - 1);
    for (int i = 0; i < ans.size(); ++i) {
        int temp = 0;
        for (int j = 0; j <= i; ++j) {
            int element_p = get(j, P);
            int element_q = get(i-j, Q);
            if (element_p != 0 && element_q != 0){
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
bch_code::summing_polynomial(const std::vector<int> &P, const std::vector<int> &Q) const {
    std::vector<int> ans(std::max(P.size(), Q.size()));
    for (int i = 0; i < ans.size(); ++i) {
        ans[i] = get(i, P) ^ get(i, Q);
    }
    return ans;
}

void bch_code::scalar_mul_vector(int scalar, std::vector<int> &vector, std::vector<int> &result) const {
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

galois_field bch_code::build_galois_field(int _n) {
    //Нужно, чтобы n|(q^m - 1)(в нашем случае q == 2, Q = q^m)
    int m;
    int max_m = 20;
    int Q;
    bool existBCH = false;

    for (m = 1, Q = 2; m <= max_m; ++m, Q <<= 1) {
        if ((Q - 1) % _n == 0) {
            existBCH = true;
            break;
        }
    }

    if (!existBCH) {
        throw std::overflow_error("BCH code doesn't exist or order is too big");
    }

//TODO Прикрутить нормальную генерацию порождающих многочленов
    if (m == 8) {
        return galois_field((1 << 8) | (1 << 4) | (1 << 3) | 5);
    }

    //Полином x^m + x^i + 1 - не всегда примитивный
    for (int i = 1; i < m; ++i) {
        int irreducible_polynomial = (1 << m) | (1 << i) | 1;
        //Нужно проверить является ли этот многочлен примитивным
        //Для этого нужно проверить, что все его корни - порождающие элементы мультипликативной группы поля
        //То есть порядок любого корня должен быть равен порядку мультпликативной группы расширенного поля, т.е. (1 << m) - 1
        try {
            return galois_field(irreducible_polynomial);
        } catch (std::invalid_argument &e) {
            continue;
        }
    }

    throw std::runtime_error("Can't create galois field");
}

[[nodiscard]] std::vector<int> bch_code::get_generating_polynomial(int _delta) const {

    std::set<std::vector<int>> minimal_polynomials;

    //если b == 1 - код БЧХ в узком смысле
    int b = 1;
    //заполнение множества минимальных многочленов
    //требуется найти минимальные многочлены для элементов alpha^b, alpha^(b+1) ... alpha^(b + delta - 2)
    for (int i = b; i <= b + _delta - 2; ++i) {
        std::vector<int> minimal_polynomial{1};
        int m_b = GF.m;
        for (int j = 1; j <= GF.m; ++j) {
            if (GF.getElement(i * (1 << j)) == GF.getElement(i)) {
                m_b = j;
                break;
            }
        }
        for (int j = 1, degree = 2; j <= m_b; ++j, degree *= 2) {
            std::vector<int> multiplier{GF.getElement(i * degree), 1};
            minimal_polynomial = multiply_polynomial(minimal_polynomial, multiplier);
        }
        minimal_polynomials.insert(minimal_polynomial);
    }

    std::vector<int> generating_polynomial{1};
    for (auto &t: minimal_polynomials) {
        generating_polynomial = multiply_polynomial(generating_polynomial, t);
    }

    return generating_polynomial;
}

bch_code::bch_code(int _n, int _delta) :
        GF(build_galois_field(_n)),
        generatingPolynomial(get_generating_polynomial(_delta)) {
    n = _n;
    k = n - generatingPolynomial.size() + 1;
    if (_delta > n) {
        throw std::runtime_error("Can't create code with delta > n");
    }
    delta = _delta;

}

std::vector<int> bch_code::shiftLeft(const std::vector<int> &a, size_t shift) {
    std::vector<int> result(a.size() + shift, 0);
    for (int i = 0; i < a.size(); ++i) {
        result[i + shift] = a[i];
    }
    return result;
}

std::vector<int> bch_code::mod(std::vector<int> a, std::vector<int> &b) const {
    while (a.size() >= b.size()) {
        std::vector<int> t = shiftLeft(b, a.size() - b.size());
        int coeff = GF.multiply_elements(a[a.size() - 1], GF.getInverse(t[t.size() - 1]));
        scalar_mul_vector(coeff, t, t);
//        t = multiply_polynomial({coeff}, t);
        t = summing_polynomial(a, t);
        shrink(t);
        a = t;
    }

    return a;
}

std::vector<int> bch_code::div(std::vector<int> a, std::vector<int> &b) const {
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

        t = multiply_polynomial({coeff}, t);
        a = summing_polynomial(t, a);
        shrink(a);
    }
    return result;
}

int bch_code::degree(std::vector<int> &polynomial) const {
    for (int i = polynomial.size() - 1; i >= 0; --i) {
        if (polynomial[i] != 0) {
            return i;
        }
    }
    return -INF;

}

std::vector<int> bch_code::gcd(std::vector<int> a, std::vector<int> b) const {
    shrink(a);
    shrink(b);

    while (!b.empty()) {
        a = mod(a, b);
        std::swap(a, b);
    }
    return a;
}

int bch_code::evaluate(const std::vector<int> &polynomial, int a) const {
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

std::vector<int> bch_code::code_word(std::vector<int> &information_word) {
    auto shifted = shiftLeft(information_word, n - k);
//    shrink(shifted);
    auto m = mod(shifted, generatingPolynomial);
    auto res = summing_polynomial(m, shifted);
    return res;
}

[[nodiscard]] std::vector<int> bch_code::get_syndrome_vector(const std::vector<int> &corrupted_word) const {
    std::vector<int> syndrome(delta - 1, 0);
    for (int i = 1; i <= 1 + delta - 2; ++i) {
        for (int j = 0; j < corrupted_word.size(); ++j) {
            if (corrupted_word[j] != 0) {
                syndrome[i - 1] ^= GF.getElement(i * j);
            }
        }
    }
    return syndrome;
}


[[nodiscard]] std::vector<int> bch_code::decoder_berlekamp_massey(const std::vector<int> &syndrome) const {
    std::vector<int> locators{1};
    std::vector<int> B{1};
    int r = 1;
    int m = 0;
    int L = 0;
    while (r <= delta - 1) {
        int delta_r = 0;
        for (int j = 0; j <= L; ++j) {
            delta_r = GF.sum_elements(
                    delta_r,
                    GF.multiply_elements(locators[j],
                                         syndrome[r - 1 - j])
            );
        }
        if (delta_r != 0) {
            std::vector<int> B_temp(shiftLeft(B, r - m));
            scalar_mul_vector(delta_r, B_temp, B_temp);
            std::vector<int> T(summing_polynomial(locators, B_temp));
            if (2 * L <= r - 1) {
                int x = GF.getInverse(delta_r);
                scalar_mul_vector(x, locators, B);
                locators = T;
                L = r - L;
                m = r;
            } else {
                locators = T;
            }
        }
        r++;
    }
    return locators;
}

std::vector<int> bch_code::get_roots(std::vector<int> &polynomial) const {
    std::vector<int> roots(polynomial.size() - 1, GF.q);
    int roots_count = 0;
    size_t polynomial_degree = polynomial.size() - 1;
    for (int i = 0; i < GF.q; ++i) {
        if (evaluate(polynomial, GF.log_to_element[i]) == 0) {
            roots[roots_count] = GF.log_to_element[i];
            roots_count++;
            if (roots_count == polynomial_degree) {
                break;
            }
        }
    }
    roots.resize(roots_count);
    return roots;
}

std::vector<int> bch_code::find_errors(std::vector<int> &corrupted_word) const {
    auto locators = decoder_berlekamp_massey(get_syndrome_vector(corrupted_word));
    auto roots = get_roots(locators);
    for (int &root: roots) {
        root = GF.getLog(GF.getInverse(root));
    }
    std::sort(roots.begin(), roots.end());
    return roots;
}

void bch_code::invert_symbol(std::vector<int> &code, int index) {
    if (code[index] == 0) {
        code[index] = 1;
    } else {
        code[index] = 0;
    }
}

std::vector<int> bch_code::fix_errors(std::vector<int> &corrupted_word) const {
    auto errors_positions = find_errors(corrupted_word);
    std::vector<int> fixed(corrupted_word);
    for (auto t: errors_positions) {
        invert_symbol(fixed, t);
    }
    return fixed;
}
