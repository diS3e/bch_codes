//
// Created by S3 on 18.05.2023.
//

#ifndef BCH_CODES_CHASE_DECODER_H
#define BCH_CODES_CHASE_DECODER_H

#include <random>
#include <chrono>
#include <algorithm>
#include "bch_codes.h"
#include "random.h"
#include<set>
#include <bitset>
#include <map>

typedef std::pair<std::pair<std::vector<int>, std::vector<int>>, std::pair<std::vector<int>, std::vector<int>>> grobner_basis;

struct chase_decoder {

    bch_code bchCode;
    struct random rnd;


    explicit chase_decoder(bch_code &bchCode1) : bchCode(bchCode1) {}

    std::vector<double> get_corrupt(std::vector<int> &random_codeword, double arg) {
        std::vector<double> ans(bchCode.n, 0);
        for (int i = 0; i < bchCode.n; ++i) {
            ans[i] = ((random_codeword[i] * 2) - 1) +
                     sqrt(static_cast<double>(bchCode.n) / (2 * arg * bchCode.k)) * rnd.normal(rnd.rng);
        }
        return ans;
    }

    static std::vector<int> get_signs(std::vector<double> &corrupted) {
        std::vector<int> ans(corrupted.size());
        for (int i = 0; i < corrupted.size(); ++i) {
            ans[i] = (corrupted[i] > 0) ? 1 : 0;
        }
        return ans;
    }

    static std::vector<double> get_reliability(std::vector<double> &corrupted) {
        std::vector<double> ans(corrupted.size());
        for (int i = 0; i < ans.size(); ++i) {
            ans[i] = std::abs(corrupted[i]);
        }
        return ans;
    }


    static std::vector<int> get_unreliable_positions(std::vector<double> &reliability, int d) {
        std::vector<std::pair<double, int>> elements;
        for (int i = 0; i < reliability.size(); ++i) {
            elements.emplace_back(reliability[i], i);
        }
        auto comparator =
                [](const std::pair<double, int> &a, const std::pair<double, int> &b) { return a.first < b.first; };
        std::sort(elements.begin(), elements.end(), comparator);
        std::vector<int> result;
        result.reserve(d);
        for (int i = 0; i < d; ++i) {
            result.push_back(elements[i].second);
        }
        return result;
    }

    static std::vector<int> intToBinary(int number) {
        std::vector<int> binary_representation;
        while (number > 0) {
            binary_representation.push_back(number % 2);
            number /= 2;
        }

        return binary_representation;
    }

    std::vector<std::vector<int>> chase_decoding_2(std::vector<double> &corrupted_word, int tau) const {
        auto signs = get_signs(corrupted_word);
        auto reliability = get_reliability(corrupted_word);
        std::vector<std::vector<int>> results;
        results.reserve(1 << 6);
        long double mu_min = std::numeric_limits<long double>::max();
        auto unreliable_positions = get_unreliable_positions(reliability, tau);
        for (int i = 0; i < (1 << tau); ++i) {
            std::vector<int> test_vector(signs);
            for (int j = 0; j < tau; ++j) {
                if (((1 << j) & i) != 0) {
                    bch_code::invert_symbol(test_vector, unreliable_positions[j]);
                }
            }
            auto fixed = bchCode.fix_errors(test_vector);

            long double mu = 0;
            for (int j = 0; j < fixed.size(); ++j) {
                mu += reliability[j];
            }
            if (mu <= mu_min) {
                if (mu < mu_min) {
                    mu_min = mu;
                    results.clear();
                }
                results.push_back(fixed);
            }
        }
        return results;
    }

    int get_even(std::vector<int> &polynomial, int degree) {
        return bchCode.get(degree * 2, polynomial);
    }

    int get_odd(std::vector<int> &polynomial, int degree) {
        return bchCode.get(degree * 2 + 1, polynomial);
    }


    std::vector<int> get_modified_syndrome(std::vector<int> &syndrome, int t) {
        std::vector<int> a(t);
        a[0] = get_even(syndrome, 0);
        for (int i = 1; i < t; ++i) {
            a[i] = get_even(syndrome, i);
            for (int j = 1; j <= i; ++j) {
                a[i] = bchCode.GF.sum_elements(a[i], bchCode.GF.multiply_elements(a[i - j], get_odd(syndrome, j - 1)));
            }
        }
        return a;
    }


    std::pair<std::pair<std::vector<int>, std::vector<int>>, std::pair<std::vector<int>, std::vector<int>>>
    fitzpatrick(std::vector<int> &syndrome, int t) {
        std::vector<std::vector<int>> b{{0},
                                        {1}};
        std::vector<int> alpha{1, -100};
        int i = 1;
        int j = 1;
        int d = 1;
        for (int k = 0; k < t; ++k) {

            alpha[j] = bch_code::get(k, bchCode.multiply_polynomial(syndrome, b[j]));
            if (alpha[i] != 0) {
                int coeff = bchCode.GF.multiply_elements(alpha[1 ^ i], bchCode.GF.getInverse(alpha[i]));
                b[i ^ 1] = bchCode.summing_polynomial(b[i ^ 1],
                                                      bchCode.multiply_polynomial({coeff}, b[i]));
                b[i] = bch_code::shiftLeft(b[i], 1);
                j = i ^ 1;
                d--;
                if (d == 0) {
                    i = i ^ 1;
                    d = 1;
                }
            } else {
                b[i ^ 1] = bch_code::shiftLeft(b[i ^ 1], 1);
                j = i;
                d++;
            }
        }
        std::vector<std::vector<int>> a(2);
        for (int k = 0; k < 2; ++k) {
            a[k] = bchCode.multiply_polynomial(syndrome, b[k]);
            a[k].resize(t);
        }
        return {{a[0], b[0]},
                {a[1], b[1]}};
    }

    std::vector<int> subs_square(std::vector<int> &polynomial) {
        std::vector<int> result(polynomial.size() << 1, 0);
        for (int i = 0; i < polynomial.size(); ++i) {
            result[i << 1] = polynomial[i];
        }
        return result;
    }

    std::vector<int> mu(std::pair<std::vector<int>, std::vector<int>> &polynomial_pair) {
        return bchCode.summing_polynomial(subs_square(polynomial_pair.second),
                                          bch_code::shiftLeft(subs_square(polynomial_pair.first), 1));
    }


    double wdeg(std::pair<std::vector<int>, std::vector<int>> &polynomial_pair, double w) {
        double deg_first = bchCode.degree(polynomial_pair.first);
        double deg_second = bchCode.degree(polynomial_pair.second);
        bool is_zero_first = deg_first == bch_code::INF;
        bool is_zero_second = deg_second == -bch_code::INF;
        if (is_zero_first && is_zero_second) {
            return 0;
        } else if (is_zero_first) {
            return deg_second + w;
        } else if (is_zero_second) {
            return deg_first;
        } else
            return std::max(deg_first, deg_second + w);
    }

    std::pair<std::pair<std::vector<int>, std::vector<int>>, std::pair<std::vector<int>, std::vector<int>>>
    kotter_iteration(double w,
                     std::pair<std::pair<std::vector<int>, std::vector<int>>, std::pair<std::vector<int>, std::vector<int>>> &G,
                     const std::vector<int> &h1, const std::vector<int> &h2, int error_locator) {
        std::vector<std::pair<std::vector<int>, std::vector<int>>> Gnew(2);
        auto g1 = G.first;
        auto g2 = G.second;
        int J = -1;
        galois_field &GF = bchCode.GF;
        std::vector<int> delta(2, 0);
        int inverse = GF.getInverse(error_locator);
        int hh1 = bchCode.evaluate(h1, inverse);
        int hh2 = bchCode.evaluate(h2, inverse);

        int inverse_square = GF.getPower(error_locator, -2);
        if (hh1 != 0) {
            int coef = GF.multiply_elements(hh2, GF.getInverse(hh1));
            delta[0] = GF.sum_elements(bchCode.evaluate(G.first.first, inverse_square),
                                       GF.multiply_elements(coef, bchCode.evaluate(G.first.second,
                                                                                   inverse_square)));
            delta[1] = GF.sum_elements(bchCode.evaluate(G.second.first, inverse_square),
                                       GF.multiply_elements(coef, bchCode.evaluate(G.second.second,
                                                                                   inverse_square)));
        } else {
            delta[0] = bchCode.evaluate(G.first.second, inverse_square);
            delta[1] = bchCode.evaluate(G.second.second, inverse_square);
        }
        std::vector<int> A;
        for (int i = 0; i < 2; ++i) {
            if (delta[i] != 0) {
                A.push_back(i);
            } else {
                Gnew[i] = (i == 0) ? G.first : G.second;
            }
        }
        if (A.size() == 1) {
            J = A[0];
        } else {
            double wdeg1 = wdeg(G.first, w);
            double wdeg2 = wdeg(G.second, w);

            if (wdeg(G.first, w) < wdeg(G.second, w)) {
                J = 0;
            } else {
                J = 1;
            }
        }
        auto &G_J = (J == 0) ? G.first : G.second;
        for (int j: A) {
            auto &G_j = (j == 0) ? G.first : G.second;
            if (j != J) {
                int coef = GF.multiply_elements(delta[J], GF.getInverse(delta[j]));
                Gnew[j].first = bchCode.summing_polynomial(bchCode.multiply_polynomial({coef}, G_j.first), G_J.first);
                Gnew[j].second = bchCode.summing_polynomial(bchCode.multiply_polynomial({coef}, G_j.second),
                                                            G_J.second);
            } else {
                Gnew[j].first = bchCode.multiply_polynomial({inverse_square, 1}, G_J.first);
                Gnew[j].second = bchCode.multiply_polynomial({inverse_square, 1}, G_J.second);
            }
        }
        return {Gnew[0], Gnew[1]};
    }

    void dfs(std::vector<int> &tree, std::vector<std::set<int>> &layers, int u, int tau, int depth) {
        layers[depth].insert(u);
        for (int i = 0; i < (tau); ++i) {
            if (((1 << i) & u) == 0 && tree[(1 << i) | u] == -1) {
                tree[(1 << i) | u] = u;
                dfs(tree, layers, (1 << i) | u, tau, depth + 1);
            }
        }
    }


    std::pair<std::vector<int>, std::vector<std::set<int>>> get_decoding_tree(int tau) {
        std::vector<int> tree(1 << tau, -1);
        std::vector<std::set<int>> layers(tau + 1, std::set<int>());


        dfs(tree, layers, 0, tau, 0);
        return {tree, layers};
    }


    std::set<std::vector<int>> fast_chase_decoding_2(std::vector<double> &corrupted_word, int tau) {
        //Инициализация
        auto signs = get_signs(corrupted_word);
        auto reliability = get_reliability(corrupted_word);

        std::set<std::vector<int>> results;

        //Поиск синдрома и модифицированного синдрома. Нахождение базиса модифицированного модуля решений
        int t = (bchCode.delta - 1) / 2;
        auto syndrome = bchCode.get_syndrome_vector(signs);
        auto modifiedSyndrome = get_modified_syndrome(syndrome, t);
        auto H = fitzpatrick(modifiedSyndrome, t);

        //Проверка решения(если ошибок <= t)

        std::vector<int> sigma;
        if (wdeg(H.first, -1) < wdeg(H.second, -1)) {
            sigma = mu(H.first);
        } else {
            sigma = mu(H.second);
        }
        auto roots = bchCode.get_roots(sigma);
        if (bchCode.degree(sigma) == roots.size()) {
            auto copy = signs;
            for(auto t : roots) {
                copy[ bchCode.GF.getLog(bchCode.GF.getInverse(t))] ^= 1;
            }
            results.insert(copy);
        }






        double w = 2 * bchCode.degree(H.second.second) - ((bchCode.delta - 1) / 2) - 0.5;
        auto caph1 = mu(H.first);
        auto caph2 = mu(H.second);

        auto unreliable_positions = get_unreliable_positions(reliability, tau);

        //Сюда идем, если не получилось продекодировать раньше. Т.е. нужно перед этим проверить правильность декодирования
        auto pair = get_decoding_tree(tau);
        auto tree = pair.first;
        auto layers = pair.second;
        std::map<int, grobner_basis> previous_layer;
        previous_layer[0] = {{{1}, {0}},
                                    {{0}, {1}}};
        for (int tree_layer = 1; tree_layer <= tau; ++tree_layer) {
            for (auto next_node: layers[tree_layer]) {
                int new_error_locator = bchCode.GF.getElement(unreliable_positions[(int) log2(next_node ^ tree[next_node])]);
                auto chase_result = kotter_iteration(w, previous_layer[tree[next_node]], caph1, caph2, new_error_locator);
                //Поиск ошибок и проверка условия декодирования
                std::vector<int> g_0;
                std::vector<int> g_1;
                int m_w = -1;
                if (wdeg(chase_result.first, w) < wdeg(chase_result.second, w)) {
                    g_0 = chase_result.first.first;
                    g_1 = chase_result.first.second;
                    m_w = 0;
                } else {
                    g_0 = chase_result.second.first;
                    g_1 = chase_result.second.second;
                    m_w = 1;
                }

                auto gcd = (bchCode.degree(g_0) > bchCode.degree(g_1)) ? bchCode.gcd(g_0, g_1) : bchCode.gcd(g_1, g_0);
                auto f_1 = bchCode.div(g_0, gcd);
                auto f_2 = bchCode.div(g_1, gcd);
                std::set<int> errors;
                for (int i = 1; i < bchCode.GF.q; ++i) {
                    int inverse = bchCode.GF.getInverse(i);
                    int inverse_square = bchCode.GF.getPower(inverse, 2);
                    int a = bchCode.evaluate(f_1, inverse_square);
                    int b = bchCode.evaluate(f_2, inverse_square);
                    int c = bchCode.evaluate(caph1, inverse);
                    int d = bchCode.evaluate(caph2, inverse);
                    if (bchCode.GF.multiply_elements(a, c) == bchCode.GF.multiply_elements(b, d)) {
                        errors.insert(i);
                    }
                }

                int delta = std::max(2 * bchCode.degree(f_1) + 2 * bchCode.degree(H.first.first) + 1, 2 * bchCode.degree(f_2) + 2 * bchCode.degree(H.second.second));

                if (delta == errors.size()) {
                    auto copy = signs;
                    for(auto t : errors) {
                        copy[bchCode.GF.getLog(t)] ^= 1;
                    }
                    results.insert(copy);
                }

            }
        }
        return results;
    }

};

#endif //BCH_CODES_CHASE_DECODER_H
