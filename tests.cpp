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
    galua_field GF(19);
    EXPECT_EQ(14, GF.getElement(11));
    EXPECT_EQ(1, GF.getElement(0));
    EXPECT_EQ(2, GF.getElement(1));
    EXPECT_EQ(4, GF.getElement(2));
    EXPECT_EQ(GF.getElement(5), GF.getElement(20));
    EXPECT_EQ(GF.getElement(5), GF.getElement(-25));
}

TEST(default_operation, getLog) {
    galua_field GF(19);
    EXPECT_EQ(11, GF.getLog(14));
    EXPECT_EQ(0, GF.getLog(1));
    EXPECT_EQ(1, GF.getLog(2));
    EXPECT_EQ(2, GF.getLog(4));
    EXPECT_EQ(5, GF.getLog(GF.getElement(5)));

}

TEST(default_operation, sum) {
    galua_field GF(19);
    EXPECT_EQ(1, GF.sum_elements(0, 1));
    EXPECT_EQ(7, GF.sum_elements(11, 12));
    EXPECT_EQ(1, GF.sum_elements(GF.getElement(8), GF.getElement(2)));

}

TEST(default_operation, multiply) {
    galua_field GF(19);
    EXPECT_EQ(3, GF.multiply_elements(3, 1));
    EXPECT_EQ(6, GF.multiply_elements(3, 2));
    EXPECT_EQ(GF.log_to_element[6], GF.multiply_elements(GF.log_to_element[10], GF.log_to_element[11]));
    EXPECT_EQ(12, GF.multiply_elements(GF.log_to_element[10], GF.log_to_element[11]));
    EXPECT_EQ(12, GF.multiply_elements(7, 14));
    EXPECT_EQ(GF.getElement(8), GF.multiply_elements(GF.getElement(5), GF.getElement(3)));
    EXPECT_EQ(0, GF.multiply_elements(0, GF.getElement(3)));
}

TEST(default_operation, simple_expression) {
    galua_field GF(19);

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

TEST(generating_polynomial, check_count) {
    auto GF = galua_field(19);
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
//    auto GF = galua_field(19);
//    std::vector<int> answer{1, 1, 1, 0, 1, 0, 0, 0, 1};
    int answer = 1 + 16 + 64 + 128 + 256;
    int generated = get_generating_polynomial(15, 7, 5);
    std::cout << std::bitset<15>(generated);
//    EXPECT_EQ(answer.size(), generated.size());
//    for (int i = 0; i < answer.size(); ++i) {
//        EXPECT_EQ(answer[i], generated[i]);
//    }
    EXPECT_EQ(answer, generated);
}

TEST(generating_polynomial, check_generating_polynomial_2) {
//    auto GF = galua_field(19);
//    std::vector<int> answer{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int answer = (1 << 15) - 1;
    int generated = get_generating_polynomial(15, 7, 9);
    std::cout << std::bitset<15>(generated);

//    EXPECT_EQ(answer.size(), generated.size());
//    for (int i = 0; i < answer.size(); ++i) {
//        EXPECT_EQ(answer[i], generated[i]);
//    }
    EXPECT_EQ(answer, generated);
}


TEST(generating_polynomial, coding_word) {
    int generated = get_generating_polynomial(15, 7, 5);
//    std::vector<int> information_word{1, 1, 0, 0, 1};
    int information_word = 25;
    int coded_word = code_word(information_word, generated, 15, 7);
    EXPECT_EQ(coded_word, (1 << 12) + (1 << 11) + (1<<8) + (1<<7) + (1<<6) +32 + 16 + 4 + 2);

}

