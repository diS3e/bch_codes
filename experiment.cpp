//
// Created by S3 on 25.05.2023.
//

#include <iostream>
#include <iomanip>
#include "chase_decoder.h"

struct random rnd;

bool words_equals(std::vector<int> const &a, std::vector<int> const &b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

void collect_data(int n, int delta, int precision, int tau) {
    bch_code bchCode(n, delta);
    chase_decoder chaseDecoder(bchCode);

    for (double Eb = 4.0; Eb < 6.0; Eb += 0.5) {
        double correct = 0;
        double all = precision;
        for (int tries = 0; tries < all; tries++) {
            auto information_word = rnd.get_random_word(bchCode.k);
            auto coded_word = bchCode.code_word(information_word);
            std::vector<double> corrupt = chaseDecoder.get_corrupt(coded_word, pow(10.0, (Eb / 10)));
            auto results = chaseDecoder.chase_decoding_2(corrupt, tau);
            for (auto &t: results) {
                if (words_equals(t, coded_word)) {
                    correct++;
                    break;
                }
            }
        }
        std::cout << std::fixed << std::setprecision(6) << Eb << ";" << correct / all << '\n';
    }
}


void collect_mean_numbers_vector(int n, int delta, int precision, int tau) {
    bch_code bchCode(n, delta);
    chase_decoder chaseDecoder(bchCode);

    for (double Eb = 0.0; Eb < 6.0; Eb += 0.5) {
        double correct = 0;
        double all = precision;
        for (int tries = 0; tries < all; tries++) {
            auto information_word = rnd.get_random_word(bchCode.k);
            auto coded_word = bchCode.code_word(information_word);
            std::vector<double> corrupt = chaseDecoder.get_corrupt(coded_word, pow(10.0, (Eb / 10)));
            auto results = chaseDecoder.chase_decoding_2(corrupt, tau);
            std::set<std::vector<int>> set(results.begin(), results.end());
            correct += set.size();

        }
        std::cout << std::fixed << std::setprecision(6) << Eb << ";" << correct / all << '\n';
    }
}


int main() {
//    collect_mean_numbers_vector(15, 5, 10000, 3);
    bch_code bchCode(255, 21);
    chase_decoder chaseDecoder(bchCode);
    int all = 1000;
    int correct = 0;
    for (int qqq = 0; qqq < all; ++qqq) {
        auto informationWord = rnd.get_random_word(bchCode.k);
        auto codedWord = bchCode.code_word(informationWord);
        std::set<int> input_errors;
        int last_error = -1;
        while (input_errors.size() != ((bchCode.delta - 1) / 2) + 1) {
            last_error = rnd.rnd(0, bchCode.n - 1);
            input_errors.insert(last_error);
        }


        auto corruptedWord = codedWord;
        for (auto t: input_errors) {
            corruptedWord[t] ^= 1;
        }

        auto syndrome = bchCode.get_syndrome_vector(corruptedWord);
        auto modifiedSyndrome = chaseDecoder.get_modified_syndrome(syndrome);
        auto H = chaseDecoder.fitzpatrick(modifiedSyndrome);
        double w = 2 * bchCode.degree(H.second.second) - ((bchCode.delta - 1) / 2) - 0.5;
        auto caph1 = chaseDecoder.mu(H.first);
        auto caph2 = chaseDecoder.mu(H.second);

        std::pair<std::pair<std::vector<int>, std::vector<int>>, std::pair<std::vector<int>, std::vector<int>>> begin_basis = {{{1}, {0}},
                                                                                                                               {{0}, {1}}};
        auto chase_result = begin_basis;
        chase_result = chaseDecoder.kotter_iteration(w, begin_basis, caph1, caph2,
                                                     bchCode.GF.getElement(rnd.rnd(0, bchCode.n - 1)));

        std::vector<int> g_0;
        std::vector<int> g_1;
        int m_w = -1;
        if (chaseDecoder.wdeg(chase_result.first, w) < chaseDecoder.wdeg(chase_result.second, w)) {
            g_0 = chase_result.first.first;
            g_1 = chase_result.first.second;
            m_w = 0;
        } else {
            g_0 = chase_result.second.first;
            g_1 = chase_result.second.second;
            m_w = 1;

        }

        auto gcd = (bchCode.degree(g_0) > bchCode.degree(g_1)) ? bchCode.gcd(g_0, g_1) : bchCode.gcd(g_1, g_0);
        auto f_1 = bchCode.div(g_0, gcd);
        auto f_2 = bchCode.div(g_1, gcd);
        std::set<int> errors;
        for (int i = 1; i < bchCode.GF.q; ++i) {
            int inverse = bchCode.GF.getInverse(i);
            int inverse_square = bchCode.GF.getPower(inverse, 2);
            int a = bchCode.evaluate(f_1, inverse_square);
            int b = bchCode.evaluate(f_2, inverse_square);
            int c = bchCode.evaluate(caph1, inverse);
            int d = bchCode.evaluate(caph2, inverse);
            if (bchCode.GF.multiply_elements(a, c) == bchCode.GF.multiply_elements(b, d)) {
                errors.insert(i);
            }
        }

        if (input_errors.size() == errors.size()) {
            bool flag = true;
            for (auto t: input_errors) {
                if (errors.count(bchCode.GF.getElement(t)) == 0) {
                    flag = false;
                }
            }
            if (flag) {
                correct++;
            }
        }
    }
    std::cout << "All: " << all << "\n";
    std::cout << "Correct: " << correct << "\n";

}
