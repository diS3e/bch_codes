#ifndef BCH_CODES_POLYNOMIALS_FIELD_H
#define BCH_CODES_POLYNOMIALS_FIELD_H

#include <vector>
#include "galois_field.h"

struct polynomials_field {
    static constexpr int INF = 100000;

    galois_field GF;

    explicit polynomials_field(galois_field galoisField) : GF(galoisField){}

    static int get(int i, const std::vector<int> &v);

    static void shrink(std::vector<int> &a);

    void scalar_mul_vector(int scalar, std::vector<int> &vector, std::vector<int> &result) const;

    static std::vector<int> shiftLeft(const std::vector<int> &a, std::size_t shift);

    [[nodiscard]] std::vector<int> mul(const std::vector<int> &P, const std::vector<int> &Q) const;

    [[nodiscard]] std::vector<int> sum(const std::vector<int> &P, const std::vector<int> &Q) const;

    static int degree(std::vector<int> &polynomial) ;

    [[nodiscard]] std::vector<int> gcd(std::vector<int> a, std::vector<int> b) const;

    std::vector<int> mod(std::vector<int> a, std::vector<int> &b) const;

    std::vector<int> div(std::vector<int> a, std::vector<int> &b) const;

    [[nodiscard]] int evaluate(const std::vector<int> &polynomial, int a) const;

    std::vector<int> get_roots(std::vector<int> &polynomial) const;
};

#endif //BCH_CODES_POLYNOMIALS_FIELD_H
