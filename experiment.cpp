#include <iostream>
#include <iomanip>
#include "chase_decoder.h"
#include <set>
#include <algorithm>

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

void printVector(const std::string &prefix, const std::vector<int> &vector) {
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
//std::vector<double> corrupt{0.666238, 0.440147, -0.246487, -0.187352, 1.8618, 1.34823, 0.790773, 0.0331272, -0.169725, -1.10385, -0.741989, 0.0566463, 1.91129, -0.48681, -0.954932};
//            printVector("Test:", corrupt);

//            std::cout << "Test" << ": [";
//            for (int i = 0; i < corrupt.size(); ++i) {
//                std::cout << corrupt[i] << ", ";
//            }
//            std::cout << "]\n";
//
            auto signs = chase_decoder::get_signs(corrupt);
            auto reliab = chase_decoder::get_reliability(corrupt);
            auto unrel = chase_decoder::get_unreliable_positions(reliab, tau);
            auto results_fast = chaseDecoder.fast_chase_decoding_2(corrupt, tau);
            std::cout << "Fast:" << results_fast.size() << std::endl;

            auto results_default = chaseDecoder.chase_decoding_2(corrupt, tau);
            std::cout << "Default:" << results_default.size() << std::endl;

            if (std::count(results_fast.begin(), results_fast.end(), coded_word) == 0 &&
                std::count(results_default.begin(), results_default.end(), coded_word) != 0) {
                printVector("\nCoded word", coded_word);
                printVector("Corrupted ", signs);
                printVector("Unreliable", unrel);
                std::cout << "Test" << ": [";
                for (int i = 0; i < corrupt.size(); ++i) {
                    std::cout << corrupt[i] << ", ";
                }
                std::cout << "]\n";
                std::cout << "Fast result:\n";
                for (auto &t: results_fast) {
                    printVector("", t);
                    printVector("Syndrome", bchCode.get_syndrome_vector(t));

                }
                std::cout << "Default:\n";
                for (auto &t: results_default) {
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
//    std::cout << "\nRuntime =  << clock()/1000.0 << "\n";


//    std::vector<int> strange_correct_vector{1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0};
//    std::vector<int> strange_vector        {1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1};
//
//    auto errors = bchCode.find_errors(strange_vector);
    std::cout << (compare_two_methods(15, 5, 6000, 1) ? "EQUAL" : "NOT EQUAL");
//    collect_data_with_fast_chase(63, 11, 10000, 6);

}
