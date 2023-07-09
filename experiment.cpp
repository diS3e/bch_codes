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

    for (double Eb = -2.0; Eb < 6.0; Eb += 0.5) {
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


void collect_data_with_fast_chase(int n, int delta, int precision, int tau) {
    bch_code bchCode(n, delta);
    chase_decoder chaseDecoder(bchCode);

    for (double Eb = -2.0; Eb < 6.0; Eb += 0.5) {
        double correct = 0;
        double all = precision;
        for (int tries = 0; tries < all; tries++) {
            auto information_word = rnd.get_random_word(bchCode.k);
            auto coded_word = bchCode.code_word(information_word);
            std::vector<double> corrupt = chaseDecoder.get_corrupt(coded_word, pow(10.0, (Eb / 10)));
            auto results = chaseDecoder.fast_chase_decoding_2(corrupt, tau);
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

void printVector(const std::string& prefix, const std::vector<int> &vector) {
    std::cout << prefix << ": [";
    for (int i = 0; i < vector.size(); ++i) {
        std::cout << vector[i] << ", ";
    }
    std::cout << "]\n";
}

bool compare_two_methods(int n, int delta, int precision, int tau) {
    bch_code bchCode(n, delta);
    chase_decoder chaseDecoder(bchCode);

    for (double Eb = 0.0; Eb < 6.0; Eb += 0.5) {
        double correct = 0;
        double all = precision;
        for (int tries = 0; tries < all; tries++) {
            auto information_word = rnd.get_random_word(bchCode.k);
            auto coded_word = bchCode.code_word(information_word);
            std::vector<double> corrupt = chaseDecoder.get_corrupt(coded_word, pow(10.0, (Eb / 10)));
            auto signs = chase_decoder::get_signs(corrupt);
            auto reliab = chase_decoder::get_reliability(corrupt);
            auto unrel = chase_decoder::get_unreliable_positions(reliab, tau);
            auto results_fast = chaseDecoder.fast_chase_decoding_2(corrupt, tau);


            auto results_default = chaseDecoder.chase_decoding_2(corrupt, tau);

            if (results_fast.count(coded_word) == 0 &&
                std::count(results_default.begin(), results_default.end(), coded_word) != 0) {
                printVector("\nCoded word", coded_word);
                printVector("Corrupted ", signs);
                printVector("Unreliable", unrel);
                std::cout << "Fast result:\n";
                for (auto &t: results_fast) {
                    printVector("", t);
                    printVector("Syndrome", bchCode.get_syndrome_vector(t));

                }
                std::cout << "Default:\n";
                for (std::vector<int> &t: results_default) {
                    printVector("", t);
                    printVector("Syndrome", bchCode.get_syndrome_vector(t));
                }
                return false;
            }
        }
        std::cout << std::fixed << std::setprecision(6) << Eb << ";" << correct / all << '\n';
    }
    return true;
}

int main() {
    bch_code bchCode(15, 5);
    chase_decoder chaseDecoder(bchCode);


//    std::cout << "Fast results:\n";
//    collect_data_with_fast_chase(63, 11, 10000, 6);
//    std::cout << "\nRuntime = " << clock()/1000.0 << "\n";
//
//    std::cout << "Default results:\n";
//    collect_data(63,11,10000, 6);
//    std::cout << "\nRuntime = " << clock()/1000.0 << "\n";


//    std::vector<int> strange_correct_vector{1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0};
//    std::vector<int> strange_vector        {1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1};
//
//    auto errors = bchCode.find_errors(strange_vector);
    std::cout << (compare_two_methods(15, 5, 1000, 1) ? "EQUAL" : "NOT EQUAL");


}
