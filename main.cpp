//TODO декодер чейза

//TODO графики такие же как для порядковых статистик

//TODO среднее количество различных векторов

#include <iostream>
#include <fstream>
//#include <iomanip>
#include "chase_decoder.h"

random rnd;

std::vector<int> strToVector(std::string &string) {
    std::vector<int> result(string.size());
    for (int i = 0; i < string.size(); i++) {
        result[i] = string[i] - '0';
    }
    return result;
}

void printVector(std::ostream &out, std::vector<int> const &vector) {
    for (auto t: vector) {
        out << t;
    }
    out << '\n';
}

bool words_equals(std::vector<int> const &a, std::vector<int> const &b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}


int main(int argc, char *argv[]) {

    if (argc == 1) {
        std::cout << "No command line arguments\n";
    }

    auto start_time = std::chrono::steady_clock::now();
    int n = std::stoi(argv[1]);
    int delta = std::stoi(argv[2]);
    bch_code bchCode(n, delta);
    std::string mode = argv[3];

    std::cout << "BchCode(" << bchCode.n << ", " << bchCode.k << ", " << bchCode.delta << ")\n\n";
    if (mode == "Files" || mode == "-f") {
        std::string operation = argv[4];
        std::ifstream in(argv[5]);
        std::ofstream out(argv[6]);
        std::string word;
        while (in >> word) {
            auto vector = strToVector(word);
            if (operation == "--encoding") {
                if (vector.size() == bchCode.k) {
                    auto encoded_word = bchCode.code_word(vector);
                    printVector(out, encoded_word);
                } else {
                    out << "Invalid size of encoding word\n";
                }
            } else if (operation == "--decoding") {
                if (vector.size() == n) {
                    auto decoded_word = bchCode.fix_errors(vector);
                    for (int i = n - bchCode.k; i < n; i++) {
                        out << decoded_word[i];
                    }
                    out << '\n';
                } else {
                    out << "Invalid size of decoding word\n";
                }
            }
        }
    } else if (mode == "Numbers" || mode == "-n") {
        int errors = std::stoi(argv[4]);
        int iterations = std::stoi(argv[5]);
        int successfully_decoded = 0;
        for (int i = 0; i < iterations; ++i) {
            auto information_word = rnd.get_random_word(bchCode.k);
            auto coded_word = bchCode.code_word(information_word);
            auto corrupted_word(coded_word);
            std::set<int> set;

            for (int j = 0; j < errors; ++j) {
                int x = rnd.rnd(0, n - 1);
                while (set.count(x) != 0) {
                    x = rnd.rnd(0, n - 1);
                }
                bch_code::invert_symbol(corrupted_word, x);
                set.insert(x);
            }

//            for (int j = 0; j < errors; ++j) {
//                if ()
//            }
            auto fixed_word = bchCode.fix_errors(corrupted_word);
            if (words_equals(fixed_word, coded_word)) {
                successfully_decoded++;
            }
        }

        std::cout << "Number of errors:\t" << errors << std::endl;
        std::cout << "Total tests:\t\t" << iterations << std::endl;
        std::cout << "Successfully decoded:\t" << successfully_decoded << std::endl;
    } else if (mode == "Probability" || mode == "-p") {
        double probability = std::stod(argv[4]);

        int iterations = std::stoi(argv[5]);
        int successfully_decoded = 0;
        for (int i = 0; i < iterations; ++i) {
            auto information_word = rnd.get_random_word(bchCode.k);
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
    return 0;
}
