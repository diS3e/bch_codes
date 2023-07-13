#ifndef BCH_CODES_BCH_CODES_H
#define BCH_CODES_BCH_CODES_H

#include "polynomials_field.h"

struct bch_code {
public:
    polynomials_field field;
    galois_field const & GF;
    std::vector<int> generatingPolynomial;
    int n;
    int k;
    int delta;

    bch_code(int _n, int _delta);

    static galois_field build_galois_field(int _n);

    std::vector<int> find_errors(std::vector<int> &corrupted_word) const;

    std::vector<int> code_word(std::vector<int> &information_word);

    std::vector<int> fix_errors(std::vector<int> &corrupted_word) const;

//private:
    [[nodiscard]] std::vector<int> get_generating_polynomial(int _delta) const;

    [[nodiscard]] std::vector<int> get_syndrome_vector(const std::vector<int> &corrupted_word) const;

    [[nodiscard]] std::vector<int> decoder_berlekamp_massey(const std::vector<int> &syndrome) const;

    friend struct chase_decoder;
};

#endif //BCH_CODES_BCH_CODES_H
