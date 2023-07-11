//
// Created by S3 on 18.05.2023.
//

#ifndef BCH_CODES_CHASE_DECODER_H
#define BCH_CODES_CHASE_DECODER_H

#include "bch_codes.h"
#include "random.h"
#include <vector>

typedef std::pair<std::pair<std::vector<int>, std::vector<int>>, std::pair<std::vector<int>, std::vector<int>>> grobner_basis;

struct chase_decoder {

    bch_code bchCode;
    struct random rnd;

    explicit chase_decoder(bch_code &bchCode1);

    std::vector<double> get_corrupt(std::vector<int> &random_codeword, double arg);

    static std::vector<int> get_signs(std::vector<double> &corrupted);

    static std::vector<double> get_reliability(std::vector<double> &corrupted);


    static std::vector<int> get_unreliable_positions(std::vector<double> &reliability, int d);

    static std::vector<int> intToBinary(int number);

    std::vector<std::vector<int>> chase_decoding_2(std::vector<double> &corrupted_word, int tau);

    int get_even(std::vector<int> &polynomial, int degree);

    int get_odd(std::vector<int> &polynomial, int degree);


    std::vector<int> get_modified_syndrome(std::vector<int> &syndrome, int t);


    grobner_basis fitzpatrick(std::vector<int> &syndrome, int t);


    double wdeg(std::pair<std::vector<int>, std::vector<int>> &polynomial_pair, double w);

    grobner_basis kotter_iteration(double w, grobner_basis &G,
                                   const std::vector<int> &h1, const std::vector<int> &h2, int error_locator);

    void dfs(std::vector<int> &tree, std::vector<std::vector<int>> &layers, int u, int tau, int depth);


    std::pair<std::vector<int>, std::vector<std::vector<int>>> get_decoding_tree(int tau);


    std::vector<std::vector<int>> fast_chase_decoding_2(std::vector<double> &corrupted_word, int tau);

};

#endif //BCH_CODES_CHASE_DECODER_H
