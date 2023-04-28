//
// Created by s3 on 23.04.23.
//

#ifndef BCH_CODES_BCH_CODES_H
#define BCH_CODES_BCH_CODES_H

#include <cmath>
#include <map>
#include <bitset>


struct galua_field {
    std::map<int, int> log_to_element;
    std::map<int, int> element_to_log;
    int m;
    int GF_power;

    explicit galua_field(int generating_polynomial) {
        int prev = 1;
        m = (int) log2(generating_polynomial);
        GF_power = (1 << m);
        for (int i = 0; i < (1 << m) - 1; ++i) {
            std::cout << "alpha^" << i << ": " << std::bitset<5>(prev) << "\n";
            log_to_element.insert({i, prev});
            element_to_log.insert({prev, i});
            prev <<= 1;
            if (prev >= (1 << m)) {
                prev = prev ^ generating_polynomial;
            }
        }
    }

    void check_containing_element(int a) {
        if (!(0 <= a && a < (1 << m))) {
            throw "Element not exists in this field";
        }
    }

//Складывает элементы поля
    int sum_elements(int a, int b) {
        check_containing_element(a);
        check_containing_element(b);
        return a ^ b;
    }

//Умножает элементы поля
    int multiply_elements(int a, int b) {
        check_containing_element(a);
        check_containing_element(b);
        if (a == 0 || b == 0) {
            return 0;
        }
        return log_to_element[(element_to_log[a] + element_to_log[b]) % (GF_power - 1)];
    }
//Возвращает логарифм элемента
    int getLog(int element) {
        check_containing_element(element);
        if (element == 0) {
            throw "Can't get log from 0";
        }
        return element_to_log[element];
    }
//Возвращает элемент по логарифму
    int getElement(int log) {
        log = ((log % (GF_power - 1)) + (GF_power - 1)) % (GF_power - 1);
        return log_to_element[log];
    }
};

int get(int i, const std::vector<int>& v){
    return (i < v.size()) ? v[i] : 0;
}

std::vector<int> multiply_polynomial(galua_field GF, const std::vector<int> &P, const std::vector<int> &Q) {
    std::vector<int> ans(P.size() + Q.size() - 1);
    for (int i = 0; i < ans.size(); ++i) {
        int temp = 0;
        for (int j = 0; j <= i; ++j) {
            temp = GF.sum_elements(temp, GF.multiply_elements(get(j, P), get(i-j, Q)));
        }
        ans[i] = temp;
    }
    return ans;
}

int get_generating_polynomial(int n, int k, int delta) {
    //Нужно, чтобы n|(q^m - 1)(в нашем случае q == 2, Q = q^m)
    int m;
    int max_m = 10;
    int Q;
    bool existBCH = false;
    for (m = 1, Q = 2; m <= max_m; ++m, Q *= 2) {
        if ((Q - 1) % n == 0) {
            existBCH = true;
            break;
        }
    }
    if (!existBCH) {
        throw "BCH code doesn't exist or order is too big";
    }

    //многочлен вида x^m + x + 1 = (10...011) неприводим над GF(2)
    galua_field GF = galua_field(1 << m | 3);

    std::set<std::vector<int>> minimal_polynomials;

    //если b == 1 - код БЧХ в узком смысле
    int b = 1;

    //заполнение множества минимальных многочленов
    //требуется найти минимальные многочлены для элементов alpha^b, alpha^(b+1) ... alpha^(b + delta - 2)
    for (int i = b; i <= b + delta - 2; ++i) {
        std::vector<int> minimal_polynomial{1};
        int m_b = m;
        for (int j = 1; j <= m; ++j) {
            if (GF.getElement(i * (1 << j)) == GF.getElement(i)) {
                m_b = j;
                break;
            }
        }
        for (int j = 1, degree = 2; j <= m_b; ++j, degree *= 2) {
            std::vector<int> multiplier{1, GF.getElement(i * degree)};
            minimal_polynomial = multiply_polynomial(GF, minimal_polynomial, multiplier);
        }
        minimal_polynomials.insert(minimal_polynomial);
    }

    std::vector generating_polynomial{1};
    for(auto &t: minimal_polynomials) {
        generating_polynomial = multiply_polynomial(GF, generating_polynomial, t);
    }

    int result_polynomial = 0;
    for(int i = 0; i < generating_polynomial.size(); i++) {
        std::cout << generating_polynomial[i] << ' ';
        result_polynomial |= (1 << i) * generating_polynomial[generating_polynomial.size() - i - 1];
    }

    return result_polynomial;
}

int code_word(int information_word, int generated_polynomial, int n, int k) {
    information_word <<= (n - k);
    std::cout << "\nmoved:\t\t\t\t" << std::bitset<15>(information_word) << '\n';
    std::cout << "\ngenerated:\t\t\t" << std::bitset<15>(generated_polynomial) << '\n';
//    information_word % generated_polynomial;
    std::cout << "moved with mod:\t\t" << std::bitset<15>(generated_polynomial ^ (information_word % generated_polynomial))<< '\n';
    std::cout << "result:\t\t\t\t" << std::bitset<15>((information_word % generated_polynomial) + information_word)<< '\n';

    return (information_word % generated_polynomial) & information_word;

}

#endif //BCH_CODES_BCH_CODES_H
