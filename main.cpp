//TODO декодер чейза
//TODO графики такие же как для порядковых статистик
//TODO среднее количество различных векторов

#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
//#include <random>
#include "bch_codes.h"
//#include "bch_codes.cpp"
#include "galois_field.h"
//#include "galois_field.cpp"
#include "random.h"
#include "chase_decoder.h"

random rnd;

std::vector<int> strToVector(std::string &string) {
    std::vector<int> result(string.size());
    for (int i = 0; i < string.size(); i++) {
        result[i] = string[i] - '0';
    }
    return result;
}

void printVector(std::ostream &out, std::vector<int> const&vector) {
    for (auto t: vector) {
        out << t;
    }
    out << '\n';
}


std::vector<int> intToBinary(int number) {
    std::vector<int> binary_representation;
    while (number > 0) {
        binary_representation.push_back(number % 2);
        number /= 2;
    }

    return binary_representation;
}

bool words_equals(std::vector<int> const &a, std::vector<int> const &b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

void collect_data() {

        std::uniform_int_distribution<int> range(0, 1);

    bch_code bchCode(63, 35, 12);
    chase_decoder chaseDecoder(bchCode);

    for (double Eb = 5.0; Eb < 6.0; Eb += 0.1) {
        double correct = 0;
        double all = 1000;
        for(int tries = 0; tries < all; tries++){
            auto information_word = rnd.get_random_word(bchCode.k);
            auto coded_word = bchCode.code_word(information_word);
            std::vector<double> corrupt = chaseDecoder.get_corrupt(coded_word, pow(10.0, (Eb / 10)));
            auto results = chaseDecoder.chase_decoding_2(corrupt, 6);
            for(auto &t: results) {
                if (words_equals(t, coded_word)) {
                    correct++;
                    break;
                }
            }
        }
        std::cout << std::fixed << std::setprecision(5) << Eb << ";" << correct / all << '\n';
    }
}

int main(int argc, char *argv[]) {

    if (argc == 1) {
        std::cout << "No command line arguments\n";
    }

    auto start_time = std::chrono::steady_clock::now();
    int n = std::stoi(argv[1]);
    int k = std::stoi(argv[2]);
    int delta = std::stoi(argv[3]);
    bch_code bchCode(n, k, delta);
    std::string mode = argv[4];

    std::cout << "BchCode(" << n << ", " << k << ", " << delta << ")\n\n";
    if (mode == "Files" || mode == "-f") {
        std::string operation = argv[5];
        std::ifstream in(argv[6]);
        std::ofstream out(argv[7]);
        std::string word;
        while (in >> word) {
            auto vector = strToVector(word);
            if (operation == "--encoding") {
                if (vector.size() == k) {
                    auto encoded_word = bchCode.code_word(vector);
                    printVector(out, encoded_word);
                } else {
                    out << "Invalid size of encoding word\n";
                }
            } else if (operation == "--decoding") {
                if (vector.size() == n) {
                    auto decoded_word = bchCode.fix_errors(vector);
                    for (int i = n - k; i < n; i++) {
                        out << decoded_word[i];
                    }
                    out << '\n';
                } else {
                    out << "Invalid size of decoding word\n";
                }
            }
        }
    } else if (mode == "Numbers" || mode == "-n") {
        int errors = std::stoi(argv[5]);
        int iterations = std::stoi(argv[6]);
        int successfully_decoded = 0;
        for (int i = 0; i < iterations; ++i) {
            auto information_word = rnd.get_random_word(k);
            auto coded_word = bchCode.code_word(information_word);
            auto corrupted_word(coded_word);
            for (int j = 0; j < errors; ++j) {
                bch_code::invert_symbol(corrupted_word, rnd.rnd(1, n - 1));
            }
            auto fixed_word = bchCode.fix_errors(corrupted_word);
            if (words_equals(fixed_word, coded_word)) {
                successfully_decoded++;
            }
        }

        std::cout << "Number of errors:\t" << errors << std::endl;
        std::cout << "Total tests:\t\t" << iterations << std::endl;
        std::cout << "Successfully decoded:\t" << successfully_decoded << std::endl;
    } else if (mode == "Probability" || mode == "-p") {
        double probability = std::stod(argv[5]);

        int iterations = std::stoi(argv[6]);
        int successfully_decoded = 0;
        for (int i = 0; i < iterations; ++i) {
            auto information_word = rnd.get_random_word(k);
            auto coded_word = bchCode.code_word(information_word);
            auto corrupted_word(coded_word);

            for (int j = 0; j < corrupted_word.size(); ++j) {
                if (rnd.uniform_real(rnd.rng) <= probability) {
                    bch_code::invert_symbol(corrupted_word, j);
                }
            }

            auto fixed_word = bchCode.fix_errors(corrupted_word);
            if (words_equals(fixed_word, coded_word)) {
                successfully_decoded++;
            }
        }

        std::cout << "Error probability:\t" << probability << std::endl;
        std::cout << "Total tests:\t\t" << iterations << std::endl;
        std::cout << "Incorrectly decoded:\t" << 1 - ((double) successfully_decoded / iterations) << std::endl;
    } else {
        std::cerr << "Unexpected mode token\n";
        return 0;
    }

    auto end_time = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "\nThe time: " << elapsed_ms.count() << "ms\n";
//    std::vector<double> reliability{1, 2, 3, -1, 5, 6};
//    printVector(std::cout, chase_decoder::get_unreliable_positions(reliability, 2));
    return 0;
//    collect_data();
}
