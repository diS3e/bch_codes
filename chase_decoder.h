//
// Created by S3 on 18.05.2023.
//

#ifndef BCH_CODES_CHASE_DECODER_H
#define BCH_CODES_CHASE_DECODER_H

#include <random>
#include <chrono>
#include <algorithm>
#include "bch_codes.h"
#include "random.h"
#include<set>

struct chase_decoder {
    bch_code bchCode;
    struct random rnd;


    explicit chase_decoder(bch_code &bchCode1) : bchCode(bchCode1) {}

    std::vector<double> get_corrupt(std::vector<int> &random_codeword, double arg) {
        std::vector<double> ans(bchCode.n, 0);
        for (int i = 0; i < bchCode.n; ++i) {
            ans[i] = ((random_codeword[i] * 2) - 1) +
                     sqrt(static_cast<double>(bchCode.n) / (2 * arg * bchCode.k)) * rnd.normal(rnd.rng);
        }
        return ans;
    }

    static std::vector<int> get_signs(std::vector<double> &corrupted) {
        std::vector<int> ans(corrupted.size());
        for (int i = 0; i < corrupted.size(); ++i) {
            ans[i] = (corrupted[i] > 0) ? 1 : 0;
        }
        return ans;
    }

    static std::vector<double> get_reliability(std::vector<double> &corrupted) {
        std::vector<double> ans(corrupted.size());
        for (int i = 0; i < ans.size(); ++i) {
            ans[i] = std::abs(corrupted[i]);
        }
        return ans;
    }


    static std::vector<int> get_unreliable_positions(std::vector<double> &reliability, int d) {
        std::vector<std::pair<double, int>> elements;
        for (int i = 0; i < reliability.size(); ++i) {
            elements.emplace_back(reliability[i], i);
        }
        auto comparator =
                [](const std::pair<double, int> &a, const std::pair<double, int> &b) { return a.first < b.first; };
        std::sort(elements.begin(), elements.end(), comparator);
        std::vector<int> result;
        result.reserve(d);
        for (int i = 0; i < d; ++i) {
            result.push_back(elements[i].second);
        }
        return result;
    }

    static std::vector<int> intToBinary(int number) {
        std::vector<int> binary_representation;
        while (number > 0) {
            binary_representation.push_back(number % 2);
            number /= 2;
        }

        return binary_representation;
    }

    std::vector<std::vector<int>> chase_decoding_2(std::vector<double> &corrupted_word, int tau) const {
        auto signs = get_signs(corrupted_word);
        auto reliability = get_reliability(corrupted_word);
        std::vector<std::vector<int>> results;
        results.reserve(1 << 6);
        long double mu_min = std::numeric_limits<long double>::max();
        auto unreliable_positions = get_unreliable_positions(reliability, tau);
        for (int i = 0; i < (1 << tau); ++i) {
            std::vector<int> test_vector(signs);
            for (int j = 0; j < tau; ++j) {
                if (((1 << j) & i) != 0) {
                    bch_code::invert_symbol(test_vector, unreliable_positions[j]);
                }
            }
            auto fixed = bchCode.fix_errors(test_vector);

            long double mu = 0;
            for (int j = 0; j < fixed.size(); ++j) {
                mu += reliability[j];
            }
            if (mu <= mu_min) {
                if (mu < mu_min) {
                    mu_min = mu;
                    results.clear();
                }
                results.push_back(fixed);
            }
        }
        return results;
    }
};

#endif //BCH_CODES_CHASE_DECODER_H
