//
// Created by S3 on 08.05.2023.
//

#include <cmath>
#include <stdexcept>
#include "galois_field.h"

galois_field::galois_field(int generating_polynomial) {
    int prev = 1;
    m = (int) log2(generating_polynomial);
    log_to_element = new int[(1 << m) - 1];
    element_to_log = new int[(1 << m) - 1];
    for (int i = 0; i < (1 << m) - 1; ++i) {
        log_to_element[i] = prev;
        element_to_log[prev] = i;
        prev <<= 1;
        if (prev >= (1 << m)) {
            prev = prev ^ generating_polynomial;
        }
        if (prev == 1 && i != (1 << m) - 2) {
            throw std::invalid_argument("Polynomial isn't primitive");
        }
    }
}

void galois_field::check_containing_element(int a) const {
    if (!(0 <= a && a < (1 << m))) {
        throw std::range_error("Element not exists in this field");
    }
}

[[nodiscard]] int galois_field::sum_elements(int a, int b) const {
    check_containing_element(a);
    check_containing_element(b);
    return a ^ b;
}

[[nodiscard]] int galois_field::multiply_elements(int a, int b) const {
    check_containing_element(a);
    check_containing_element(b);
    if (a == 0 || b == 0) {
        return 0;
    }
    return log_to_element[(element_to_log[a] + element_to_log[b]) % ((1 << m) - 1)];
}

[[nodiscard]] int galois_field::getLog(int element) const {
    check_containing_element(element);
    if (element == 0) {
        throw std::range_error("Can't get log from 0");
    }
    return element_to_log[element];
}

[[nodiscard]] int galois_field::getElement(int log) const {
    log = ((log % ((1 << m) - 1)) + ((1 << m) - 1)) % ((1 << m) - 1);
    return log_to_element[log];
}

[[nodiscard]] int galois_field::getInverse(int a) const {
    return getElement(-getLog(a));
}

[[nodiscard]] int galois_field::getPower(int base, int pow) const {
    if (base == 0) {
        if (pow == 0) {
            return 1;
        } else {
            return 0;
        }
    } else {
        return getElement(getLog(base) * pow);
    }
}