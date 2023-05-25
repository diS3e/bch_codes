//
// Created by s3 on 23.04.23.
//

#ifndef BCH_CODES_BCH_CODES_H
#define BCH_CODES_BCH_CODES_H

#include <vector>
#include "galois_field.h"


struct bch_code {
public:
    galois_field GF;
    std::vector<int> generatingPolynomial;
    int n;
    int k;
    int delta;
    bch_code(int _n, int _delta);
    static int get(int i, const std::vector<int> &v);
    static void invert_symbol(std::vector<int>& code, int index) ;
    static galois_field build_galois_field(int _n) ;
    static std::vector<int> shrink(const std::vector<int> &a);
    static std::vector<int> shiftLeft(const std::vector<int> &a, std::size_t shift) ;

    std::vector<int> find_errors(std::vector<int> &corrupted_word) const ;

    std::vector<int> code_word(std::vector<int> &information_word) ;

    std::vector<int> fix_errors(std::vector<int> &corrupted_word) const ;

private:
    [[nodiscard]] std::vector<int> multiply_polynomial(const std::vector<int> &P, const std::vector<int> &Q) const;

    [[nodiscard]] std::vector<int> summing_polynomial(const std::vector<int> &P, const std::vector<int> &Q) const ;

    [[nodiscard]] std::vector<int> get_generating_polynomial(int _delta) const ;

    std::vector<int> mod(std::vector<int> a, std::vector<int> &b) const ;

    [[nodiscard]] int evaluate(std::vector<int> polynomial, int a) const ;

    [[nodiscard]] std::vector<int> get_syndrome_vector(const std::vector<int> &corrupted_word) const ;

    [[nodiscard]] std::vector<int> decoder_berlekamp_massey(const std::vector<int> &syndrome) const ;

    std::vector<int> get_roots(std::vector<int> &polynomial) const ;
};

#endif //BCH_CODES_BCH_CODES_H
