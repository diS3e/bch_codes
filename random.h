//
// Created by S3 on 18.05.2023.
//

#ifndef BCH_CODES_RANDOM_H
#define BCH_CODES_RANDOM_H

#include <random>
#include <chrono>

struct random {
    std::mt19937 rng;
    std::uniform_int_distribution<int> uniform_int;
    std::uniform_real_distribution<> uniform_real;
    std::normal_distribution<double> normal;



    random() : rng((std::chrono::steady_clock::now().time_since_epoch().count())),
               uniform_int(0, std::numeric_limits<int>::max()),
               uniform_real(0.0, 1.0),
               normal(0, 1) {

    }

    int rnd(int l, int r) {
        return (uniform_int(rng) % (r - l + 1) + l);
    }

    std::vector<int> get_random_word(int length) {
        std::vector<int> result(length);
        for (auto &t: result) t = rnd(0, 1);
            return result;
    }


};

#endif //BCH_CODES_RANDOM_H
