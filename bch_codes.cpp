//
// Created by s3 on 23.04.23.
//
#include "bch_codes.h"
#include <algorithm>
#include <set>

int bch_code::get(int i, const std::vector<int> &v) {
    return (0 <= i && i < v.size()) ? v[i] : 0;
}

// Метод shrink - убирается нулевые старшие степени
std::vector<int> bch_code::shrink(const std::vector<int> &a) {
    std::vector<int> result(a);
    while (!result.empty() && result[result.size() - 1] == 0) {
        result.erase(result.end() - 1);
    }
    return result;
}

[[nodiscard]] std::vector<int> bch_code::multiply_polynomial(const std::vector<int> &P, const std::vector<int> &Q) const {
    std::vector<int> ans(P.size() + Q.size() - 1);
    for (int i = 0; i < ans.size(); ++i) {
        int temp = 0;
        for (int j = 0; j <= i; ++j) {
            temp = GF.sum_elements(temp, GF.multiply_elements(get(j, P), get(i - j, Q)));
        }
        ans[i] = temp;
    }
    return ans;
}

[[nodiscard]] std::vector<int> bch_code::summing_polynomial(const std::vector<int> &P, const std::vector<int> &Q) const {
    std::vector<int> ans(std::max(P.size(), Q.size()));
    for (int i = 0; i < ans.size(); ++i) {
        ans[i] = GF.sum_elements(get(i, P), get(i, Q));
    }
    return ans;
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

    //Полином x^m + x^i + 1
    for (int i = 1; i < m; ++i) {
        int irreducible_polynomial = (1 << m) | (1 << i) | 1;
        //Нужно проверить является ли этот многочлен примитивным
        //Для этого нужно проверить, что все его корни - порождающие элементы мультипликативной группы поля
        //То есть порядок любого корня должен быть равен порядку мультпликативной группы расширенного поля, т.е. (1 << m) - 1
        try {
            return  galois_field(irreducible_polynomial);
        } catch (std::invalid_argument& e) {
            continue;
        }
    }

    //многочлен вида x^m + x + 1 = (10...011) неприводим над GF(2)
//    return galois_field(1 << m | 3);
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

    std::vector generating_polynomial{1};
    for (auto &t: minimal_polynomials) {
        generating_polynomial = multiply_polynomial(generating_polynomial, t);
    }

    return generating_polynomial;
}

bch_code::bch_code(int _n, int _k, int _delta) :
        GF(build_galois_field(_n)),
        generatingPolynomial(get_generating_polynomial(_delta)) {
    n = _n;
    k = _k;
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
        a = shrink(summing_polynomial(a, t));
    }

    return a;
}

[[nodiscard]] int bch_code::evaluate(std::vector<int> polynomial, int a) const {
    int value = 0;
    for (int j = 0; j < polynomial.size(); ++j) {
        value = GF.sum_elements(value, GF.multiply_elements(polynomial[j], GF.getPower(a, j)));
    }
    return value;
}

std::vector<int> bch_code::code_word(std::vector<int> &information_word) {
    auto shifted = shiftLeft(information_word, n - k);
    auto m = mod(shifted, generatingPolynomial);
    auto res = summing_polynomial(m, shifted);
    return res;
}

[[nodiscard]] std::vector<int> bch_code::get_syndrome_vector(const std::vector<int> &corrupted_word) const {
    std::vector<int> syndrome(delta - 1);
    for (int i = 1; i <= 1 + delta - 2; ++i) {
        syndrome[i - 1] = evaluate(corrupted_word, GF.getElement(i));
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
            std::vector<int> T;
            std::vector<int> poly{delta_r};
            poly = shiftLeft(poly, r - m);
            T = summing_polynomial(locators,
                                   multiply_polynomial(poly, B));
            if (2 * L <= r - 1) {
                B = multiply_polynomial({GF.getInverse(delta_r)},
                                        locators);
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
    std::vector<int> roots;
    for (int i = 0; i < (1 << GF.m); ++i) {
        if (evaluate(polynomial, i) == 0) {
            roots.push_back(i);
            if (roots.size() == polynomial.size() - 1) {
                break;
            }
        }
    }
    return roots;
}

std::vector<int> bch_code::find_errors(std::vector<int> &corrupted_word) const {
    auto locators = decoder_berlekamp_massey(get_syndrome_vector(corrupted_word));
    auto roots = get_roots(locators);
    for (int & root : roots) {
        root = GF.getLog(GF.getInverse(root));
    }
    std::sort(roots.begin(), roots.end());
    return roots;
}

void bch_code::invert_symbol(std::vector<int>& code, int index) {
    if (code[index] == 0) {
        code[index] = 1;
    } else {
        code[index] = 0;
    }
}

std::vector<int> bch_code::fix_errors(std::vector<int> &corrupted_word) const {
    auto errors_positions = find_errors(corrupted_word);
    std::vector<int> fixed(corrupted_word);
    for(auto t: errors_positions) {
        invert_symbol(fixed, t);
    }
    return fixed;
}