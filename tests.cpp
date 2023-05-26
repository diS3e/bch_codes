//
// Created by S3 on 23.04.2023.
//
#include <gtest/gtest.h>
#include <random>
#include "bch_codes.h"
#include "chase_decoder.h"

//Нас интересуют только GF(2^m)
//Элементы этого поля задаются, как вектора над полем GF(2), т.е. двоичные представления int
//Пример:
//Примитивный многочлен: x^4 + x + 1 <-> (10011)
//Элементы
//(0000) = 0
//(0001) = 1
//(0010) = alpha
//(0100) = alpha^2
//(1000) = alpha^3
//(0011) = alpha^4 = (10000) + (10011)
//(0110) = alpha^5
//(1100) = alpha^6
//(1011) = alpha^7 = (11000) + (10011)
//(0101) = alpha^8 = (10110) + (10011)
//(1010) = alpha^9
//(0111) = alpha^10 = (10100) + (10011)
//(1110) = alpha^11
//(1111) = alpha^12 = (11100) + (10011)
//(1101) = alpha^13 = (11110) + (10011)
//(1001) = alpha^14 = (11010) + (10011)
//(0001) = alpha^15 = (10010) + (10011)

//Генерация случайных чисел
//std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
//int rnd(int l, int r) {
//    std::uniform_int_distribution<int> range(l, r);
//    return range(rng);
//}
//
//std::random_device rd;
//std::mt19937 e2(rd());


void vector_equals(std::vector<int> const &a, std::vector<int> const &b) {
    EXPECT_EQ(a.size(), b.size());
    for (int i = 0; i < a.size(); ++i) {
        EXPECT_EQ(a[i], b[i]);
    }
}

std::vector<int> intToBinary(int number) {
    std::vector<int> binary_representation;
    while (number > 0) {
        binary_representation.push_back(number % 2);
        number /= 2;
    }

    return binary_representation;
}

int rnd(int l, int r) {
    static std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<int> range(l, r);
    return range(rng);
}

void test_bch_code(bch_code& bchCode, int errors, int iterations) {
    int random_word = rnd(0, (1 << bchCode.k) - 1);
    auto information_word = intToBinary(random_word);
    information_word.resize(bchCode.k, 0);
    auto coded_word = bchCode.code_word(information_word);
    for (int j = 0; j < iterations; ++j) {
        auto copy(coded_word);
        for (int i = 0; i < errors; ++i) {
            int x = rnd(0, bchCode.n - 1);
            bch_code::invert_symbol(copy, x);
        }
        auto fixed_word = bchCode.fix_errors(copy);
        vector_equals(fixed_word, coded_word);
    }
}



void printVector(std::vector<int> &vector) {
    for (auto t: vector) {
        std::cout << t;
    }
    std::cout << '\n';
}

TEST(galois_field, getElement) {
    galois_field GF(19);
    EXPECT_EQ(14, GF.getElement(11));
    EXPECT_EQ(1, GF.getElement(0));
    EXPECT_EQ(2, GF.getElement(1));
    EXPECT_EQ(4, GF.getElement(2));
    EXPECT_EQ(GF.getElement(5), GF.getElement(20));
    EXPECT_EQ(GF.getElement(5), GF.getElement(-25));
}

TEST(galois_field, getLog) {
    galois_field GF(19);
    EXPECT_EQ(11, GF.getLog(14));
    EXPECT_EQ(0, GF.getLog(1));
    EXPECT_EQ(1, GF.getLog(2));
    EXPECT_EQ(2, GF.getLog(4));
    EXPECT_EQ(5, GF.getLog(GF.getElement(5)));

}

TEST(galois_field, sum) {
    galois_field GF(19);
    EXPECT_EQ(1, GF.sum_elements(0, 1));
    EXPECT_EQ(7, GF.sum_elements(11, 12));
    EXPECT_EQ(1, GF.sum_elements(GF.getElement(8), GF.getElement(2)));

}

TEST(galois_field, multiply) {
    galois_field GF(19);
    EXPECT_EQ(3, GF.multiply_elements(3, 1));
    EXPECT_EQ(6, GF.multiply_elements(3, 2));
    EXPECT_EQ(GF.log_to_element[6], GF.multiply_elements(GF.log_to_element[10], GF.log_to_element[11]));
    EXPECT_EQ(12, GF.multiply_elements(GF.log_to_element[10], GF.log_to_element[11]));
    EXPECT_EQ(12, GF.multiply_elements(7, 14));
    EXPECT_EQ(GF.getElement(8), GF.multiply_elements(GF.getElement(5), GF.getElement(3)));
    EXPECT_EQ(0, GF.multiply_elements(0, GF.getElement(3)));
}

TEST(galois_field, simple_expression) {
    galois_field GF(19);

    EXPECT_EQ(GF.getElement(6), GF.multiply_elements(
            GF.sum_elements(
                    GF.getElement(2),
                    GF.multiply_elements(
                            GF.getElement(20),
                            GF.getElement(3)
                    )
            ),
            GF.getElement(-9)
    ));
}

TEST(galois_field, creating_from_non_primitive) {
    EXPECT_ANY_THROW(galois_field((1 << 5) - 1));
}

TEST(galois_field, generating_in_bch_codes) {
    EXPECT_NO_THROW(bch_code bchCode1(63,  2));
    EXPECT_NO_THROW(bch_code bchCode2(63,  5));
    EXPECT_NO_THROW(bch_code bchCode3(63,  6));
    EXPECT_NO_THROW(bch_code bchCode4(63,  11));
    EXPECT_NO_THROW(bch_code bchCode5(127, 3));
}

TEST(galois_field, generating_in_bch_codes_2) {
//    EXPECT_NO_THROW(bch_code bchCode1(23, 7));
}

TEST(polynomial_operation, intToBinary) {
    int x = 19;
    auto result = intToBinary(x);
    EXPECT_EQ(result[0], 1);
    EXPECT_EQ(result[1], 1);
    EXPECT_EQ(result[2], 0);
    EXPECT_EQ(result[3], 0);
    EXPECT_EQ(result[4], 1);
}


TEST(polynomial_operation, shift) {
    int x = 19;
    auto result = bch_code::shiftLeft(intToBinary(x), 2);
    EXPECT_EQ(result[0], 0);
    EXPECT_EQ(result[1], 0);
    EXPECT_EQ(result[2], 1);
    EXPECT_EQ(result[3], 1);
    EXPECT_EQ(result[4], 0);
    EXPECT_EQ(result[5], 0);
    EXPECT_EQ(result[6], 1);
}

TEST(polynomial_operation, shrink) {
    int x = 19;
    auto result = intToBinary(x);
    result.push_back(0);
    result.push_back(0);
    result.push_back(0);
    result.push_back(0);
    EXPECT_EQ(result[0], 1);
    EXPECT_EQ(result[1], 1);
    EXPECT_EQ(result[2], 0);
    EXPECT_EQ(result[3], 0);
    EXPECT_EQ(result[4], 1);
    EXPECT_EQ(result[5], 0);
    EXPECT_EQ(result[6], 0);
    EXPECT_EQ(result[7], 0);
    result = bch_code::shrink(result);
    EXPECT_EQ(result.size(), 5);
    EXPECT_EQ(result[0], 1);
    EXPECT_EQ(result[1], 1);
    EXPECT_EQ(result[2], 0);
    EXPECT_EQ(result[3], 0);
    EXPECT_EQ(result[4], 1);
}

TEST(generating_polynomial, check_count) {
    auto GF = galois_field(19);
    EXPECT_EQ(0, GF.sum_elements(GF.getElement(4),
                                 GF.sum_elements(
                                         GF.getElement(2),
                                         GF.sum_elements(
                                                 GF.getElement(8),
                                                 GF.getElement(16)
                                         )
                                 )));
    EXPECT_EQ(0, GF.sum_elements(GF.getElement(24),
                                 GF.sum_elements(GF.getElement(6),
                                                 GF.sum_elements(
                                                         GF.getElement(12),
                                                         GF.sum_elements(GF.getElement(20),
                                                                         GF.sum_elements(GF.getElement(10),
                                                                                         GF.getElement(18))
                                                         )
                                                 )
                                 )));
    EXPECT_EQ(1, GF.sum_elements(GF.getElement(28),
                                 GF.sum_elements(
                                         GF.getElement(26),
                                         GF.sum_elements(
                                                 GF.getElement(14),
                                                 GF.getElement(22)
                                         )
                                 )));
    EXPECT_EQ(1, GF.getElement(30));
    EXPECT_EQ(1, GF.getElement(90));
}

TEST(decoding, fix_0_error) {
    bch_code bchCode(15, 5);
    test_bch_code(bchCode, 0, 10000);
}

TEST(decoding, fix_1_error) {
    bch_code bchCode(15, 5);
    test_bch_code(bchCode, 1, 10000);
}

TEST(decoding, fix_2_errors) {
    bch_code bchCode(15, 5);
    test_bch_code(bchCode, 2, 10000);
}

TEST(decoding, fix_3_errors) {
//    bch_code bchCode(23,  7);
//    test_bch_code(bchCode, 3, 1000);
}

TEST(decoding, fix_5_errors) {
    bch_code bchCode(63, 11);
    test_bch_code(bchCode, 5, 10000);
}

TEST(decoding, big_code) {
    bch_code bchCode(255, 5);
    test_bch_code(bchCode, 2, 1000);
}

bool words_equals(std::vector<int> const &a, std::vector<int> const &b) {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

TEST(chase_decoding, base) {
    bch_code bchCode(63, 11);
    chase_decoder chaseDecoder(bchCode);

    for (double Eb = 0.0; Eb < 6.0; Eb += 0.1) {
        double correct = 0;
        double all = 10;
        for(int tries = 0; tries < all; tries++){
            std::vector<int> information_word(bchCode.k);
            for (int & j : information_word) {
                j = rnd(0, 1);
            }
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
