#include <map>
#include <algorithm>
#include "chase_decoder.h"

chase_decoder::chase_decoder(bch_code &bchCode1) : bchCode(bchCode1), field(bchCode.field), GF(field.GF) {}

std::vector<double> chase_decoder::get_corrupt(std::vector<int> &random_codeword, double arg) {
    std::vector<double> ans(bchCode.n, 0);
    for (int i = 0; i < bchCode.n; ++i) {
        ans[i] = ((random_codeword[i] * 2) - 1) +
                 sqrt(static_cast<double>(bchCode.n) / (2 * arg * bchCode.k)) * rnd.normal(rnd.rng);
    }
    return ans;
}

std::vector<int> chase_decoder::get_signs(std::vector<double> &corrupted) {
    std::vector<int> ans(corrupted.size());
    for (int i = 0; i < corrupted.size(); ++i) {
        ans[i] = (corrupted[i] > 0) ? 1 : 0;
    }
    return ans;
}

std::vector<double> chase_decoder::get_reliability(std::vector<double> &corrupted) {
    std::vector<double> ans(corrupted.size());
    for (int i = 0; i < ans.size(); ++i) {
        ans[i] = std::abs(corrupted[i]);
    }
    return ans;
}


std::vector<int> chase_decoder::get_unreliable_positions(std::vector<double> &reliability, int d) {
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

std::set<std::vector<int>> chase_decoder::chase_decoding_2(std::vector<double> &corrupted_word, int tau) const {
    auto signs = get_signs(corrupted_word);
    auto reliability = get_reliability(corrupted_word);
    std::set<std::vector<int>> results;
//    results.reserve(1 << 6);
    long double mu_min = std::numeric_limits<long double>::max();
    auto unreliable_positions = get_unreliable_positions(reliability, tau);
    for (int i = 0; i < (1 << tau); ++i) {
        std::vector<int> test_vector(signs);
        for (int j = 0; j < tau; ++j) {
            if (((1 << j) & i) != 0) {
                test_vector[unreliable_positions[j]] ^= 1;
            }
        }
        auto fixed = bchCode.fix_errors(test_vector);
        auto checking_syndrome = bchCode.get_syndrome_vector(fixed);
        polynomials_field::shrink(checking_syndrome);
        long double mu = 0;
        if (checking_syndrome.empty()) {

            for (int j = 0; j < fixed.size(); ++j) {
                if (fixed[j] != signs[j]) {
                    mu += reliability[j];
                }
            }
            if (mu <= mu_min) {
                if (mu < mu_min) {
                    mu_min = mu;
                    results.clear();
                }
                results.insert(fixed);
            }
        }
    }
    return results;
}

int chase_decoder::get_even(std::vector<int> &polynomial, int degree) {
    return polynomials_field::get(degree * 2, polynomial);
}

int chase_decoder::get_odd(std::vector<int> &polynomial, int degree) {
    return polynomials_field::get(degree * 2 + 1, polynomial);
}


std::vector<int> chase_decoder::get_modified_syndrome(std::vector<int> &syndrome, int t) const {
    std::vector<int> a(t);
    a[0] = get_even(syndrome, 0);
    for (int i = 1; i < t; ++i) {
        a[i] = get_even(syndrome, i);
        for (int j = 1; j <= i; ++j) {
            a[i] ^= bchCode.field.GF.multiply_elements(a[i - j], get_odd(syndrome, j - 1));
        }
    }
    return a;
}


grobner_basis chase_decoder::fitzpatrick(std::vector<int> &syndrome, int t) const {
    std::vector<std::vector<int>> b{{0},
                                    {1}};
    std::vector<int> alpha{1, -100};
    int i = 1;
    int j = 1;
    int d = 1;
    for (int k = 0; k < t; ++k) {

        alpha[j] = polynomials_field::get(k, field.mul(syndrome, b[j]));
        if (alpha[i] != 0) {
            int coeff = bchCode.GF.multiply_elements(alpha[1 ^ i], bchCode.GF.getInverse(alpha[i]));
            b[i ^ 1] = field.sum(b[i ^ 1],
                                 field.mul({coeff}, b[i]));
            b[i] = polynomials_field::shiftLeft(b[i], 1);
            j = i ^ 1;
            d--;
            if (d == 0) {
                i = i ^ 1;
                d = 1;
            }
        } else {
            b[i ^ 1] = polynomials_field::shiftLeft(b[i ^ 1], 1);
            j = i;
            d++;
        }
    }
    std::vector<std::vector<int>> a(2);
    for (int k = 0; k < 2; ++k) {
        a[k] = field.mul(syndrome, b[k]);
        a[k].resize(t);
    }
    return {{a[0], b[0]},
            {a[1], b[1]}};
}

std::vector<int> mu(std::pair<std::vector<int>, std::vector<int>> &polynomial_pair) {
    std::vector<int> result(std::max(polynomial_pair.second.size() << 1, 1 + (polynomial_pair.first.size() << 1)));
    for (int i = 0; i < result.size(); ++i) {
        if ((i & 1) == 0) {
            result[i] = polynomials_field::get(i >> 1, polynomial_pair.second);
        } else {
            result[i] = polynomials_field::get((i - 1) >> 1, polynomial_pair.first);
        }
    }
    return result;
}


double chase_decoder::wdeg(std::pair<std::vector<int>, std::vector<int>> &polynomial_pair, double w) {
    double deg_first = polynomials_field::degree(polynomial_pair.first);
    double deg_second = polynomials_field::degree(polynomial_pair.second);
    bool is_zero_first = deg_first == -polynomials_field::INF;
    bool is_zero_second = deg_second == -polynomials_field::INF;
    if (is_zero_first && is_zero_second) {
        return 0;
    } else if (is_zero_first) {
        return deg_second + w;
    } else if (is_zero_second) {
        return deg_first;
    } else {
        return std::max(deg_first, deg_second + w);
    }
}

grobner_basis chase_decoder::kotter_iteration(double w, grobner_basis &G,
                                              const std::vector<int> &h1, const std::vector<int> &h2,
                                              int error_locator) const {
    grobner_basis Gnew;
    int J = -1;
    std::vector<int> delta(2, 0);
    int inverse = GF.getInverse(error_locator);
    int hh1 = field.evaluate(h1, inverse);
    int hh2 = field.evaluate(h2, inverse);

    int inverse_square = GF.multiply_elements(inverse, inverse);
    if (hh1 != 0) {
        int coef = GF.multiply_elements(hh2, GF.getInverse(hh1));

        delta[0] = field.evaluate(G.first.first, inverse_square) ^
                   GF.multiply_elements(coef, field.evaluate(G.first.second,
                                                             inverse_square));
        delta[1] = field.evaluate(G.second.first, inverse_square) ^
                   GF.multiply_elements(coef, field.evaluate(G.second.second,
                                                             inverse_square));
    } else {
        delta[0] = field.evaluate(G.first.second, inverse_square);
        delta[1] = field.evaluate(G.second.second, inverse_square);
    }
    std::vector<int> A;
    for (int i = 0; i < 2; ++i) {
        if (delta[i] != 0) {
            A.push_back(i);
        } else {
            if (i == 0) {
                Gnew.first = G.first;
            } else {
                Gnew.second = G.second;
            }
        }
    }
    if (A.size() == 1) {
        J = A[0];
    } else {
        if (wdeg(G.first, w) < wdeg(G.second, w)) {
            J = 0;
        } else {
            J = 1;
        }
    }
    auto &G_J = (J == 0) ? G.first : G.second;
    for (int j: A) {
        auto &G_j = (j == 0) ? G.first : G.second;
        auto &to = (j == 0) ? Gnew.first : Gnew.second;
        if (j != J) {
            int coef = GF.multiply_elements(delta[J], GF.getInverse(delta[j]));
            to.first = field.sum(field.mul({coef}, G_j.first), G_J.first);
            to.second = field.sum(field.mul({coef}, G_j.second),
                                  G_J.second);
        } else {
            to.first = field.mul({inverse_square, 1}, G_J.first);
            to.second = field.mul({inverse_square, 1}, G_J.second);
        }
    }
    return Gnew;
}

void
chase_decoder::dfs(std::vector<int> &tree, std::vector<std::vector<int>> &layers, int u, int tau, int depth) const {
    layers[depth].push_back(u);
    for (int i = 0; i < (tau); ++i) {
        if (((1 << i) & u) == 0 && tree[(1 << i) | u] == -1) {
            tree[(1 << i) | u] = u;
            dfs(tree, layers, (1 << i) | u, tau, depth + 1);
        }
    }
}


std::pair<std::vector<int>, std::vector<std::vector<int>>> chase_decoder::get_decoding_tree(int tau) const {
    std::vector<int> tree(1 << tau, -1);
    std::vector<std::vector<int>> layers(tau + 1, std::vector<int>());
    dfs(tree, layers, 0, tau, 0);
    return {tree, layers};
}

std::set<std::vector<int>> chase_decoder::fast_chase_decoding_2(std::vector<double> &corrupted_word, int tau) const {
    //Инициализация
    long double mu_min = std::numeric_limits<long double>::max();
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
    auto roots = field.get_roots(sigma);
    if (polynomials_field::degree(sigma) == roots.size()) {
        auto copy = signs;
        double mu = 0;
        for (auto root: roots) {
            mu += reliability[GF.getLog(GF.getInverse(root))];
            copy[GF.getLog(GF.getInverse(root))] ^= 1;
        }
        if (mu <= mu_min) {
            if (mu < mu_min) {
                mu_min = mu;
            }
            results.insert(copy);
        }
        results.insert(copy);
    }


    double w = 2 * polynomials_field::degree(H.second.second) - ((bchCode.delta - 1) / 2) - 0.5;
    auto caph1 = mu(H.first);
    auto caph2 = mu(H.second);

    std::vector<int> evaluated_caph1(bchCode.GF.q - 1);
    std::vector<int> evaluated_caph2(bchCode.GF.q - 1);
    for (int i = 1; i < bchCode.GF.q; ++i) {
        int inverse = bchCode.GF.getInverse(i);
        int c = field.evaluate(caph1, inverse);
        int d = field.evaluate(caph2, inverse);
        evaluated_caph1[i] = c;
        evaluated_caph2[i] = d;
    }

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
            int new_error_locator = bchCode.GF.getElement(
                    unreliable_positions[(int) log2(next_node ^ tree[next_node])]);
            auto chase_result = kotter_iteration(w, previous_layer[tree[next_node]], caph1, caph2,
                                                 new_error_locator);
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

            auto gcd = (polynomials_field::degree(g_0) > polynomials_field::degree(g_1))
                       ? field.gcd(g_0, g_1) : field.gcd(g_1, g_0);
            auto f_1 = field.div(g_0, gcd);
            auto f_2 = field.div(g_1, gcd);
            std::vector<int> errors;
            for (int i = 1; i < bchCode.GF.q; ++i) {
                int inverse = bchCode.GF.getInverse(i);
                int inverse_square = bchCode.GF.multiply_elements(inverse, inverse);
                int a = field.evaluate(f_1, inverse_square);
                int b = field.evaluate(f_2, inverse_square);
                if (bchCode.GF.multiply_elements(a, evaluated_caph1[i]) ==
                    bchCode.GF.multiply_elements(b, evaluated_caph2[i])) {
                    errors.push_back(i);
                }
            }

            int delta = std::max(2 * polynomials_field::degree(f_1) + 2 * polynomials_field::degree(H.first.first) + 1,
                                 2 * polynomials_field::degree(f_2) + 2 * polynomials_field::degree(H.second.second));

            if (delta == errors.size()) {
//                double mu = 0;
                auto copy = signs;
                for (auto item: errors) {
//                    mu += reliability[GF.getLog(item)];
                    copy[bchCode.GF.getLog(item)] ^= 1;
                }

                long double mu = 0;
                for (int j = 0; j < copy.size(); ++j) {
                    if (copy[j] != signs[j]) {
                        mu += reliability[j];
                    }
                }
                if (mu <= mu_min) {
                    if (mu < mu_min) {
                        mu_min = mu;
                        results.clear();
                    }
                    results.insert(copy);
                }
//                results.insert(copy);
            }

        }
    }
    return results;
}