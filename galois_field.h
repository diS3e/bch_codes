#ifndef BCH_CODES_GALOIS_FIELD_H
#define BCH_CODES_GALOIS_FIELD_H

struct galois_field {

    int m;
    int q;
    int *log_to_element;
    int *element_to_log;

    explicit galois_field(int generating_polynomial);


//Складывает элементы поля
    [[nodiscard]] int sum_elements(int a, int b) const;

//Умножает элементы поля
    [[nodiscard]] int multiply_elements(int a, int b) const;

//Возвращает логарифм элемента
    [[nodiscard]] int getLog(int element) const;

//Возвращает элемент по логарифму
    [[nodiscard]] int getElement(int log) const;

//Возвращает a^-1
    [[nodiscard]] int getInverse(int a) const;

//Возвращает base^pow(0^0 = 1, 0^x = 0)
    [[nodiscard]] int getPower(int base, int pow) const;
private:
    void check_containing_element(int a) const;
};

#endif //BCH_CODES_GALOIS_FIELD_H
