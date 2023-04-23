#include <iostream>
#include <bitset>
#include <cmath>
#include <map>

const int p = 2;


using namespace std;


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



map<int, int> index_to_element;
map<int, int> element_to_index;
int m;
int GF_power;


void initialize_tables(int generating_polynomial) {
    int prev = 1;
    m = (int) log2(generating_polynomial);
    GF_power = (1 << m);
    for (int i = 0; i < (1 << m) - 1; ++i) {
        cout << "alpha^" << i << ": " << bitset<4>(prev) << "\n";
        index_to_element.insert({i, prev});
        element_to_index.insert({prev, i});
        prev <<= 1;
        if (prev >= (1 << m)) {
            prev = prev ^ generating_polynomial;
        }
    }
}

//Складывает элементы поля
int sum_elements(int a, int b) {
    return a ^ b;
}

//Умножает элементы поля
int multiply_elements(int a, int b) {
    return index_to_element[(element_to_index[a] + element_to_index[b]) % (GF_power - 1)];
}

int main() {
    initialize_tables(19);
//    for(auto &t : index_to_element) {
//        cout << t.first << ' ' << t.second << '\n';
//    }

    int a = index_to_element[3];
    int b = index_to_element[2];
    cout << "result"
    return 0;
}