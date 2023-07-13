#ifndef BCH_CODES_CHASE_DECODER_H
#define BCH_CODES_CHASE_DECODER_H

#include "bch_codes.h"
#include "random.h"
#include <vector>
#include <set>
typedef std::pair<std::pair<std::vector<int>, std::vector<int>>, std::pair<std::vector<int>, std::vector<int>>> grobner_basis;

struct chase_decoder {

    bch_code bchCode;
    polynomials_field const &field;
    galois_field const &GF;
    struct random rnd;

    explicit chase_decoder(bch_code &bchCode1);

    std::vector<double> get_corrupt(std::vector<int> &random_codeword, double arg);

    static std::vector<int> get_signs(std::vector<double> &corrupted);

    static std::vector<double> get_reliability(std::vector<double> &corrupted);

    static std::vector<int> get_unreliable_positions(std::vector<double> &reliability, int d);

    std::set<std::vector<int>> chase_decoding_2(std::vector<double> &corrupted_word, int tau) const;

    std::set<std::vector<int>> fast_chase_decoding_2(std::vector<double> &corrupted_word, int tau) const;

private:
    static int get_odd(std::vector<int> &polynomial, int degree);

    static int get_even(std::vector<int> &polynomial, int degree);

    /**
     * Finding modified syndrome (11), page 7.
     * */
    std::vector<int> get_modified_syndrome(std::vector<int> &syndrome, int t) const;

    /**
     * Finding a Grobner basis for N.
     * Algorithm F, page 20.
     * */
    grobner_basis fitzpatrick(std::vector<int> &syndrome, int t) const;

    /**
     * Finding weighted degree of polynomial pair, page 3.
     * */
    static double wdeg(std::pair<std::vector<int>, std::vector<int>> &polynomial_pair, double w);

    /**
     * Finding new Grobner basis with new error locator.
     * Algorithm A. Kotter iteration on an edge of tree, page 13.
     * */
    grobner_basis kotter_iteration(double w, grobner_basis &G,
                                   const std::vector<int> &h1, const std::vector<int> &h2, int error_locator) const;

    void dfs(std::vector<int> &tree, std::vector<std::vector<int>> &layers, int u, int tau, int depth) const;

    /**
     * Creating decoding tree, page 11.
     * First is tree, second is layers array
     * */
    std::pair<std::vector<int>, std::vector<std::vector<int>>> get_decoding_tree(int tau) const;
};

#endif //BCH_CODES_CHASE_DECODER_H
