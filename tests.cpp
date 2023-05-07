//
// Created by S3 on 23.04.2023.
//
#include <gtest/gtest.h>
#include "bch_codes.h"


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

TEST(default_operation, getElement) {
    galois_field GF(19);
    EXPECT_EQ(14, GF.getElement(11));
    EXPECT_EQ(1, GF.getElement(0));
    EXPECT_EQ(2, GF.getElement(1));
    EXPECT_EQ(4, GF.getElement(2));
    EXPECT_EQ(GF.getElement(5), GF.getElement(20));
    EXPECT_EQ(GF.getElement(5), GF.getElement(-25));
}

TEST(default_operation, getLog) {
    galois_field GF(19);
    EXPECT_EQ(11, GF.getLog(14));
    EXPECT_EQ(0, GF.getLog(1));
    EXPECT_EQ(1, GF.getLog(2));
    EXPECT_EQ(2, GF.getLog(4));
    EXPECT_EQ(5, GF.getLog(GF.getElement(5)));

}

TEST(default_operation, sum) {
    galois_field GF(19);
    EXPECT_EQ(1, GF.sum_elements(0, 1));
    EXPECT_EQ(7, GF.sum_elements(11, 12));
    EXPECT_EQ(1, GF.sum_elements(GF.getElement(8), GF.getElement(2)));

}

TEST(default_operation, multiply) {
    galois_field GF(19);
    EXPECT_EQ(3, GF.multiply_elements(3, 1));
    EXPECT_EQ(6, GF.multiply_elements(3, 2));
    EXPECT_EQ(GF.log_to_element[6], GF.multiply_elements(GF.log_to_element[10], GF.log_to_element[11]));
    EXPECT_EQ(12, GF.multiply_elements(GF.log_to_element[10], GF.log_to_element[11]));
    EXPECT_EQ(12, GF.multiply_elements(7, 14));
    EXPECT_EQ(GF.getElement(8), GF.multiply_elements(GF.getElement(5), GF.getElement(3)));
    EXPECT_EQ(0, GF.multiply_elements(0, GF.getElement(3)));
}

TEST(default_operation, simple_expression) {
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
    auto result = shiftLeft(intToBinary(x), 2);
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
    result = shrink(result);
    EXPECT_EQ(result.size(), 5);
    EXPECT_EQ(result[0], 1);
    EXPECT_EQ(result[1], 1);
    EXPECT_EQ(result[2], 0);
    EXPECT_EQ(result[3], 0);
    EXPECT_EQ(result[4], 1);
}


TEST(polynomial_operation, multiply_polynomial) {
    auto GF = galois_field(19);


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


TEST(generating_polynomial, check_generating_polynomial) {
//    auto GF = galois_field(19);
    std::vector<int> answer{1, 1, 1, 0, 1, 0, 0, 0, 1};
    std::reverse(answer.begin(), answer.end());
//    int answer = 1 + 16 + 64 + 128 + 256;
    std::vector<int> generated = get_generating_polynomial(15, 7, 5).first;
//    std::cout << std::bitset<15>(generated);

    EXPECT_EQ(answer.size(), generated.size());
    for (int i = 0; i < answer.size(); ++i) {
        EXPECT_EQ(answer[i], generated[i]);
    }
//    EXPECT_EQ(answer, generated);
}

TEST(generating_polynomial, check_generating_polynomial_2) {
//    auto GF = galois_field(19);
    std::vector<int> answer{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
//    int answer = (1 << 15) - 1;
    std::vector<int> generated = get_generating_polynomial(15, 7, 9).first;
//    std::cout << std::bitset<15>(generated);

    EXPECT_EQ(answer.size(), generated.size());
    for (int i = 0; i < answer.size(); ++i) {
        EXPECT_EQ(answer[i], generated[i]);
    }
//    EXPECT_EQ(answer, generated);
}


TEST(coding, coding_word) {
    auto generated_pair = get_generating_polynomial(15, 7, 5);
    std::vector<int> generated = generated_pair.first;
//    std::vector<int> information_word{1, 1, 1, 1, 0, 0, 1};
    auto information_word = intToBinary(25 + 32 + 64);
    auto coded_word = code_word(generated_pair.second, information_word, generated, 15, 7);
    auto expected_result = intToBinary((1 << 14) + (1 << 13) + (1 << 12) +
                                       (1 << 11) + (1 << 8) + 64 + 32 + 8 + 2);
//    printVector(expected_result);
//    printVector(coded_word);
    EXPECT_EQ(expected_result.size(), coded_word.size());
    for (int i = 0; i < expected_result.size(); ++i) {
        EXPECT_EQ(expected_result[i], coded_word[i]);
    }
}

TEST(decoding, syndrome_vector) {
    auto generating_and_field = get_generating_polynomial(15, 7, 5);
    std::vector<int> generated = generating_and_field.first;
    galois_field GF = generating_and_field.second;
    auto information_word = intToBinary(25 + 32 + 64);
    auto information_word_2 = intToBinary(24 + 32 + 64);
    auto coded_word = code_word(GF, information_word, generated, 15, 7);
    auto coded_word_2 = code_word(GF, information_word, generated, 15, 7);
    std::cout << "Coded word:\n";
    printVector(coded_word);
    if (coded_word[0] == 0) {
        coded_word[0] = 1;
    } else {
        coded_word[0] = 0;
    }
    std::cout << "Corrupted_word:\n";
    printVector(coded_word);
    auto syndrome = get_syndrome_vector(GF, coded_word, 5);
    std::cout << "Syndrome:\n";
    for (int i: syndrome) {
        std::cout << syndrome[i] << ' ';
//        EXPECT_EQ(i, 0);
    }
}
void corrupt_symbol(std::vector<int>& code, int index) {
    if (code[index] == 0) {
        code[index] = 1;
    } else {
        code[index] = 0;
    }
}

TEST(decoding, decoding_word) {
    auto generating_and_field = get_generating_polynomial(15, 7, 5);
    std::vector<int> generated = generating_and_field.first;
    galois_field GF = generating_and_field.second;
    auto information_word = intToBinary(25 + 32 + 64);
    auto coded_word = code_word(GF, information_word, generated, 15, 7);

    corrupt_symbol(coded_word, 6);
    corrupt_symbol(coded_word, 3);


    auto syndrome = get_syndrome_vector(GF, coded_word, 5);
    auto locators = decoder_berlekamp_massey(GF, syndrome, 5);

    auto roots = get_roots(GF, locators);
    for (int & root : roots) {
        root = GF.getLog(GF.getInverse(root));
    }
    std::sort(roots.begin(), roots.end());
    EXPECT_EQ(roots.size(), 2);
    EXPECT_EQ(3, roots[0]);
    EXPECT_EQ(6, roots[1]);
}

//TEST(bch_code, bch_code) {
//    bch_code bchCode(15, 7, 5);
//}