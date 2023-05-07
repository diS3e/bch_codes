//
// Created by s3 on 23.04.23.
//

#ifndef BCH_CODES_BCH_CODES_H
#define BCH_CODES_BCH_CODES_H

#include <cmath>
#include <map>
#include <bitset>


struct galois_field {
    int *log_to_element;
    int *element_to_log;
    int m;

    explicit galois_field(int generating_polynomial) {
        int prev = 1;
        m = (int) log2(generating_polynomial);
        log_to_element = new int[(1 << m) - 1];
        element_to_log = new int[(1 << m) - 1];
        for (int i = 0; i < (1 << m) - 1; ++i) {
            log_to_element[i] = prev;
            element_to_log[prev] = i;
            prev <<= 1;
            if (prev >= (1 << m)) {
                prev = prev ^ generating_polynomial;
            }
        }
    }

//    galois_field(galois_field const &GF) {
//        log_to_element = new int[(1 << GF.m) - 1];
//        element_to_log = new int[(1 << GF.m) - 1];
//        memcpy(log_to_element, GF.log_to_element, (1 << GF.m) - 1);
//        memcpy(element_to_log, GF.element_to_log, (1 << GF.m) - 1);
//        m = GF.m;
//    }

//    galois_field& operator=(const galois_field &other) {
//        using std::swap;
//        if (&other != this) {
//            galois_field GF(other);
//            swap(log_to_element, GF.log_to_element);
//            swap(element_to_log, GF.element_to_log);
//            swap(m, GF.m);
//        }
//        return *this;
//    }

    void check_containing_element(int a) const {
        if (!(0 <= a && a < (1 << m))) {
            throw std::range_error("Element not exists in this field");
        }
    }

//Складывает элементы поля
    [[nodiscard]] int sum_elements(int a, int b) const {
        check_containing_element(a);
        check_containing_element(b);
        return a ^ b;
    }

//Умножает элементы поля
    [[nodiscard]] int multiply_elements(int a, int b) const {
        check_containing_element(a);
        check_containing_element(b);
        if (a == 0 || b == 0) {
            return 0;
        }
        return log_to_element[(element_to_log[a] + element_to_log[b]) % ((1 << m) - 1)];
    }

//Возвращает логарифм элемента
    [[nodiscard]] int getLog(int element) const {
        check_containing_element(element);
        if (element == 0) {
            throw std::range_error("Can't get log from 0");
        }
        return element_to_log[element];
    }

//Возвращает элемент по логарифму
    [[nodiscard]] int getElement(int log) const {
        log = ((log % ((1 << m) - 1)) + ((1 << m) - 1)) % ((1 << m) - 1);
        return log_to_element[log];
    }

    [[nodiscard]] int getInverse(int a) const {
        return getElement(-getLog(a));
    }

    [[nodiscard]] int getPower(int base, int pow) const {
        if (base == 0) {
            if (pow == 0) {
                return 1;
            } else {
                return 0;
            }
        } else {
            return getElement(getLog(base) * pow);
        }
    }
};

//Многочлен над конечным полем: vector где степени записаны по возрастанию a[i] - коэффициент при x^i
// Все операции возвращают массивы, с потенциально нулевыми старшими коэффициентами.


int get(int i, const std::vector<int> &v) {
    return (0 <= i && i < v.size()) ? v[i] : 0;
}

// Метод shrink - убирается нулевые старшие степени
std::vector<int> shrink(const std::vector<int> &a) {
    std::vector<int> result(a);
    while (!result.empty() && result[result.size() - 1] == 0) {
        result.erase(result.end() - 1);
    }
    return result;
}

std::vector<int> multiply_polynomial(const galois_field &GF, const std::vector<int> &P, const std::vector<int> &Q) {
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

std::vector<int> summing_polynomial(const galois_field &GF, const std::vector<int> &P, const std::vector<int> &Q) {
    std::vector<int> ans(std::max(P.size(), Q.size()));
    for (int i = 0; i < ans.size(); ++i) {
        ans[i] = GF.sum_elements(get(i, P), get(i, Q));
    }
    return ans;
}

std::pair<std::vector<int>, const galois_field> get_generating_polynomial(int n, int k, int delta) {
    //Нужно, чтобы n|(q^m - 1)(в нашем случае q == 2, Q = q^m)
    int m;
    int max_m = 10;
    int Q;
    bool existBCH = false;
    for (m = 1, Q = 2; m <= max_m; ++m, Q *= 2) {
        if ((Q - 1) % n == 0) {
            existBCH = true;
            break;
        }
    }
    if (!existBCH) {
        throw std::overflow_error("BCH code doesn't exist or order is too big");
    }

    //многочлен вида x^m + x + 1 = (10...011) неприводим над GF(2)
    auto GF = galois_field(1 << m | 3);

    std::set<std::vector<int>> minimal_polynomials;

    //если b == 1 - код БЧХ в узком смысле
    int b = 1;
    //заполнение множества минимальных многочленов
    //требуется найти минимальные многочлены для элементов alpha^b, alpha^(b+1) ... alpha^(b + delta - 2)
    for (int i = b; i <= b + delta - 2; ++i) {
        std::vector<int> minimal_polynomial{1};
        int m_b = m;
        for (int j = 1; j <= m; ++j) {
            if (GF.getElement(i * (1 << j)) == GF.getElement(i)) {
                m_b = j;
                break;
            }
        }
        for (int j = 1, degree = 2; j <= m_b; ++j, degree *= 2) {
            std::vector<int> multiplier{GF.getElement(i * degree), 1};
            minimal_polynomial = multiply_polynomial(GF, minimal_polynomial, multiplier);
        }
        minimal_polynomials.insert(minimal_polynomial);
    }

    std::vector generating_polynomial{1};
    for (auto &t: minimal_polynomials) {
        generating_polynomial = multiply_polynomial(GF, generating_polynomial, t);
    }

    return {generating_polynomial, GF};
}

std::vector<int> intToBinary(int number) {
    std::vector<int> binary_representation;
    while (number > 0) {
        binary_representation.push_back(number % 2);
        number /= 2;
    }

    return binary_representation;
}

std::vector<int> shiftLeft(const std::vector<int> &a, size_t shift) {
    std::vector<int> result(a.size() + shift, 0);
    for (int i = 0; i < a.size(); ++i) {
        result[i + shift] = a[i];
    }
    return result;
}

std::vector<int> mod(const galois_field &GF, std::vector<int> a, std::vector<int> &b) {
    while (a.size() >= b.size()) {
        std::vector<int> t = shiftLeft(b, a.size() - b.size());
        a = shrink(summing_polynomial(GF, a, t));
    }

    return a;
}

void printVector(std::vector<int> &vector) {
    for (auto t: vector) {
        std::cout << t;
    }
    std::cout << '\n';
}

std::vector<int>
code_word(const galois_field &GF, std::vector<int> &information_word, std::vector<int> &generated_polynomial, int n,
          int k) {
    auto shifted = shiftLeft(information_word, n - k);
    auto m = mod(GF, shifted, generated_polynomial);
    auto res = summing_polynomial(GF, m, shifted);
    return res;
}

int evaluate(const galois_field &GF, std::vector<int> polynomial, int a) {
    int value = 0;
    for (int j = 0; j < polynomial.size(); ++j) {
        value = GF.sum_elements(value, GF.multiply_elements(polynomial[j], GF.getPower(a, j)));
    }
    return value;
}

std::vector<int> get_roots(const galois_field &GF, std::vector<int> &polynomial) {
    std::vector<int> roots;
    for (int i = 0; i < (1 << GF.m); ++i) {
        if (evaluate(GF, polynomial, i) == 0) {
            roots.push_back(i);
            if (roots.size() == polynomial.size() - 1) {
                break;
            }
        }
    }
    return roots;
}

std::vector<int> get_syndrome_vector(const galois_field &GF, std::vector<int> &corrupted_word, int delta) {
    std::vector<int> syndrome(delta - 1);
    for (int i = 1; i <= 1 + delta - 2; ++i) {
        syndrome[i - 1] = evaluate(GF, corrupted_word, GF.getElement(i));
    }
    return syndrome;
}

std::vector<int> decoder_berlekamp_massey(const galois_field &GF, std::vector<int> &syndrome, int delta) {
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
            T = summing_polynomial(GF, locators,
                                   multiply_polynomial(GF, poly, B));
            if (2 * L <= r - 1) {
                B = multiply_polynomial(GF,
                                        {GF.getInverse(delta_r)},
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


//struct bch_code {
//    galois_field GF;
//    std::vector<int> generatingPolynomial;
//    int n;
//    int k;
//    int delta;
//
//
//    explicit bch_code(const std::pair<std::vector<int>, galois_field> &pair) : GF(pair.second),
//                                                                               generatingPolynomial(pair.first) {};
//
//    bch_code(int _n, int _k, int _delta) : bch_code(get_generating_polynomial(_n, _k, _delta)) {
//        n = _n;
//        k = _k;
//        delta = _delta;
//
//    };
//
//    std::vector<int> code_word(std::vector<int> &information_word) {
//        auto shifted = shiftLeft(information_word, n - k);
//        auto m = mod(GF, shifted, generatingPolynomial);
//        auto res = summing_polynomial(GF, m, shifted);
//        return res;
//    }
//
//    std::vector<int> get_syndrome_vector(const std::vector<int> &corrupted_word) {
//        std::vector<int> syndrome(delta - 1);
//        for (int i = 1; i <= 1 + delta - 2; ++i) {
//            syndrome[i - 1] = evaluate(GF, corrupted_word, GF.getElement(i));
//        }
//        return syndrome;
//    }
//
//    std::vector<int> decoder_berlekamp_massey(const std::vector<int> &syndrome) {
//        std::vector<int> locators{1};
//        std::vector<int> B{1};
//        int r = 1;
//        int m = 0;
//        int L = 0;
//        while (r <= delta - 1) {
//            int delta_r = 0;
//            for (int j = 0; j <= L; ++j) {
//                delta_r = GF.sum_elements(
//                        delta_r,
//                        GF.multiply_elements(locators[j],
//                                             syndrome[r - 1 - j])
//                );
//            }
//            if (delta_r != 0) {
//                std::vector<int> T;
//                std::vector<int> poly{delta_r};
//                poly = shiftLeft(poly, r - m);
//                T = summing_polynomial(GF, locators,
//                                       multiply_polynomial(GF, poly, B));
//                if (2 * L <= r - 1) {
//                    B = multiply_polynomial(GF,
//                                            {GF.getInverse(delta_r)},
//                                            locators);
//                    locators = T;
//                    L = r - L;
//                    m = r;
//                } else {
//                    locators = T;
//                }
//            }
//            r++;
//        }
//        return locators;
//    }
//
//    std::vector<int> find_errors(std::vector<int> &corrupted_word) {
//        auto locators = decoder_berlekamp_massey(get_syndrome_vector(corrupted_word));
//
//    }
//
//};

#endif //BCH_CODES_BCH_CODES_H