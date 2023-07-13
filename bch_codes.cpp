#include "bch_codes.h"
#include <algorithm>
#include <set>
#include <stdexcept>

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
            minimal_polynomial = field.mul(minimal_polynomial, multiplier);
        }
        minimal_polynomials.insert(minimal_polynomial);
    }

    std::vector<int> generating_polynomial{1};
    for (auto &t: minimal_polynomials) {
        generating_polynomial = field.mul(generating_polynomial, t);
    }

    return generating_polynomial;
}

bch_code::bch_code(int _n, int _delta) :
        field(build_galois_field(_n)),
        GF(field.GF),
        generatingPolynomial(get_generating_polynomial(_delta)) {
    n = _n;
    k = n - generatingPolynomial.size() + 1;
    if (_delta > n) {
        throw std::runtime_error("Can't create code with delta > n");
    }
    delta = _delta;

}

std::vector<int> bch_code::code_word(std::vector<int> &information_word) {
    auto shifted = polynomials_field::shiftLeft(information_word, n - k);
    auto m = field.mod(shifted, generatingPolynomial);
    auto res = field.sum(m, shifted);
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
            std::vector<int> B_temp(polynomials_field::shiftLeft(B, r - m));
            field.scalar_mul_vector(delta_r, B_temp, B_temp);
            std::vector<int> T(field.sum(locators, B_temp));
            if (2 * L <= r - 1) {
                int x = GF.getInverse(delta_r);
                field.scalar_mul_vector(x, locators, B);
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

std::vector<int> bch_code::find_errors(std::vector<int> &corrupted_word) const {
    auto locators = decoder_berlekamp_massey(get_syndrome_vector(corrupted_word));
    auto roots = field.get_roots(locators);
    for (int &root: roots) {
        root = GF.getLog(GF.getInverse(root));
    }
    std::sort(roots.begin(), roots.end());
    return roots;
}

std::vector<int> bch_code::fix_errors(std::vector<int> &corrupted_word) const {
    auto errors_positions = find_errors(corrupted_word);
    std::vector<int> fixed(corrupted_word);
    for (auto t: errors_positions) {
        fixed[t] ^= 1;
    }
    return fixed;
}
