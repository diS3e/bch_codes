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
        for(int tries = 0; tries < all; tries++) {
            auto information_word = rnd.get_random_word(bchCode.k);
            auto coded_word = bchCode.code_word(information_word);
            std::vector<double> corrupt = chaseDecoder.get_corrupt(coded_word, pow(10.0, (Eb / 10)));
            auto results = chaseDecoder.chase_decoding_2(corrupt, tau);
            for(auto &t: results) {
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
        for(int tries = 0; tries < all; tries++) {
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
    collect_mean_numbers_vector(63, 11, 10000, 8);
}
