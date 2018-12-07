//
// Created by Rahul Govind on 05/12/18.
//
#include "russians.h"
#include <bits/stdc++.h>

using namespace std;

void print_vector(vector<int> v) {
    for (int i = 0; i < v.size(); i++) {
        cout << v[i] << "\t";
    }
    cout << "\n";
}

void print_matrix(vector<vector<int> > v) {
    for (int i = 0; i < v.size(); i++) {
        for (int j = 0; j < v[0].size(); j++) {
            cout << v[i][j] << "\t";
        }
        cout << "\n";
    }
}

void fill_partial_edit_distance(const vector<int> &b, const vector<int> &c,
                                const vector<int> &v, const vector<int> &w,
                                vector<int> &d, vector<int> &e) {
    /*
     * b.size() == v.size() - 1
     * c.size() == w.size() - 1
     * d.size() == b.size()
     * e.size() == c.size()
     */
    unsigned long m = v.size(), n = w.size();
    vector<vector<int> > dp(m, vector<int>(n));

    // Initialize left-border and top-border
    dp[0][0] = 0;

    for (int i = 1; i < m; i++) {
        dp[i][0] = dp[i - 1][0] + b[i - 1];
    }

    for (int j = 1; j < n; j++) {
        dp[0][j] = dp[0][j - 1] + c[j - 1];
    }

    // Fill interior of dp tale
    for (int i = 1; i < m; i++) {
        for (int j = 1; j < n; j++) {
            if (v[i] == w[j]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = min(dp[i - 1][j], min(dp[i][j - 1], dp[i - 1][j - 1])) + 1;
            }
        }
    }

    for (int i = 0; i < d.size(); i++) {
        d[i] = dp[i + 1][n - 1];
    }

    for (int j = 0; j < e.size(); j++) {
        e[j] = dp[m - 1][j + 1];
    }
}

/*
 * Here, there is no extra character before v or w
 */
void fill_partial_edit_distance2(const vector<int> &b, const vector<int> &c,
                                 const vector<int> &v, const vector<int> &w,
                                 vector<int> &d, vector<int> &e) {
    /*
     * b.size() == v.size()
     * c.size() == w.size()
     * d.size() == b.size()
     * e.size() == c.size()
     */
    unsigned long m = v.size() + 1, n = w.size() + 1;
    vector<vector<int> > dp(m, vector<int>(n));

    // Initialize left-border and top-border
    dp[0][0] = 0;

    for (int i = 1; i < m; i++) {
        dp[i][0] = dp[i - 1][0] + b[i - 1];
    }

    for (int j = 1; j < n; j++) {
        dp[0][j] = dp[0][j - 1] + c[j - 1];
    }

    // Fill interior of dp tale
    for (int i = 1; i < m; i++) {
        for (int j = 1; j < n; j++) {
            if (v[i - 1] == w[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = min(dp[i - 1][j], min(dp[i][j - 1], dp[i - 1][j - 1])) + 1;
            }
        }
    }

    for (int i = 0; i < d.size(); i++) {
        d[i] = dp[i + 1][n - 1];
    }

    for (int j = 0; j < e.size(); j++) {
        e[j] = dp[m - 1][j + 1];
    }
}

vector<int> edit_distance_rightmost_col(const vector<int> &v, const vector<int> &w) {
    unsigned long m = v.size(), n = w.size();
    vector<vector<int> > dp(m, vector<int>(n));
    vector<int> result(m);

    // Initialize borders
    for (int i = 0; i < m; i++) {
        dp[i][0] = i;
    }

    for (int i = 0; i < n; i++) {
        dp[0][i] = i;
    }

    for (int i = 1; i < m; i++) {
        for (int j = 1; j < n; j++) {
            if (v[i] == w[j]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = min(dp[i - 1][j], min(dp[i][j - 1], dp[i - 1][j - 1])) + 1;
            }
        }
    }

    for (int i = 0; i < m; i++) {
        result[i] = dp[i][n - 1];
    }
    return result;
}


void fill_base_values(vector<int> &v, ulong x, int base, int offset) {
    int i = v.size() - 1;
    while (i >= 0) {
        v[i] = int(x % base) + offset;
        x /= base;
        i -= 1;
    }
}

ulong int_pow(int x, int m) {
    ulong result = 1;
    for (int i = 0; i < m; i++) {
        result *= x;
    }
    return result;
}

//void load_to_cache_map(int t, int s) {
//    auto cache_key = make_pair(t, s);
//    int *cache = cache_values(t, s);
//    ::_cache_map[cache_key] = cache;
//}

int *cache_values(int t, int s) {
    ulong tpow = int_pow(s, t);
    ulong spow = int_pow(3, t);

    ulong max_value = int_pow(s, t);
    ulong max_value2 = int_pow(3, t);

    ulong total_size = 2 * tpow * spow * max_value * max_value2 * t;
    int *cache = new int[total_size];

    cerr << "Total size: " << total_size << endl;
    int *cache_iter = cache;
    vector<int> v(t), w(t), b(t), c(t), d(t), e(t);

    for (ulong x = 0; x < max_value; x++) {
        cerr << "Completed |   x: " << x + 1 << "/" << max_value << endl;
        fill_base_values(v, x, s, 0);
        for (ulong y = 0; y < max_value; y++) {
            fill_base_values(w, y, s, 0);

            for (ulong bx = 0; bx < max_value2; bx++) {
                fill_base_values(b, bx, 3, -1);

                for (ulong cx = 0; cx < max_value2; cx++) {
                    fill_base_values(c, cx, 3, -1);

                    fill_partial_edit_distance2(b, c, v, w, d, e);

                    for (int i = 0; i < t; i++) {
                        *cache_iter = d[i];
                        cache_iter++;
                    }
                    for (int i = 0; i < t; i++) {
                        *cache_iter = e[i];
                        cache_iter++;
                    }
                }
            }
        }
    }

    return cache;
}

ulong vector_to_index(vector<int>::iterator it, vector<int>::iterator end, int base, int offset) {
    ulong result = 0;
    while (it != end) {
        result = result * base + ((*it) - offset);
        it++;
    }
    return result;
}

inline ulong vector_to_index(const int *arr, int n, int base, int offset) {
    ulong result = 0, i = 0;
    while (i < n) {
        result = result * base + (arr[i] - offset);
        i += 1;
    }
    return result;
}

vector<int> _russians(vector<int> v, vector<int> w, int s, int t, const int *cache) {
    ulong tpow = int_pow(s, t);
    ulong spow = int_pow(3, t);

    assert(v.size() % t == 1);
    assert(w.size() % t == 1);

    ulong m = v.size();
    ulong n = w.size();

    int b[16];
    int c[16];
//    vector<int> b(t), c(t);
//    vector<int> d(t), e(t);
    // Initialize
    for (int i = 0; i < t; i++) {
        b[i] = 1;
        c[i] = 1;
    }

    vector<int> rightmost_col(v.size());
    vector<int> row_above(w.size()), row_below(w.size());

    for (int j = 0; j < n; j++) {
        row_above[j] = j;
    }
    rightmost_col[0] = w.size() - 1;
    int *cache2;

    for (int i = 0; i + t < m; i += t) {
        // Starting new row. Moving down a block. Update b
        for (int k = 0; k < t; k++) {
            b[k] = 1;
        }
        ulong base_index = vector_to_index(&v[i + 1], t, s, 0);
        ulong index;
        for (int j = 0; j + t < n; j += t) {
            // Moving right. Update c
            for (int k = 0; k < t; k++) {
                c[k] = row_above[j + k + 1] - row_above[j + k];
            }

            // Get index into cache
            index = base_index;
            index = index * tpow + vector_to_index(&w[j + 1], t, s, 0);
            index = index * spow + vector_to_index(&b[0], t, 3, -1);
            index = index * spow + vector_to_index(&c[0], t, 3, -1);
            index = index * 2 * t;

            for (int k = 0; k < t; k++) {
                if (k == 0) {
                    b[0] = row_above[j] + cache[index] - row_above[j + t];
                } else {
                    b[k] = cache[index + k] - cache[index + k - 1];
                }
                if (j + 2 * t >= n) {
                    rightmost_col[i + k + 1] = row_above[n - 1 - t] + cache[index + k];
                }
            }

            index += t;
            for (int k = 0; k < t; k++) {
                row_below[j + k + 1] = row_above[j] + cache[index + k];
            }
        }

        // Starting new row after this. So row below will be the row above
        row_above = row_below;
        row_above[0] = i + t;
    }
    return rightmost_col;
}


vector<int> russians(vector<int> v, vector<int> w, int s, int t, const int *cache) {
    if (v.size() >= t + 1 && w.size() >= t + 1) {
        // We hack our way through to get russians to work for different sizes
        ulong initial_v_size = v.size();
        ulong w2_size;
        if (w.size() % t == 0) {
            w2_size = w.size() - t + 1;
        } else {
            w2_size = w.size() - w.size() % t + 1;
        }
        while (v.size() % t != 1) {
            v.push_back(0);
        }

        vector<int> w2(w.begin(), w.begin() + w2_size);
        vector<int> result_temp = _russians(v, w2, s, t, cache);
        vector<int> result = result_temp;

        // Extend this edit distance
        for (int j=w2.size(); j < w.size(); j++) {
            result[0] = j;
            for (int i=1; i < v.size(); i++) {
                if (v[i] == w[j]) {
                    result[i] = result_temp[i-1];
                } else {
                    result[i] = min(result_temp[i], min(result_temp[i-1], result[i-1])) + 1;
                }
            }
            result_temp = result;
        }
        return vector(result.begin(), result.begin() + initial_v_size);
    } else {
        cout << "V: "; print_vector(v);
        cout << "W: "; print_vector(w);
        return edit_distance_rightmost_col(v, w);
    }
//    if (v.size() >= t + 1 && w.size() >= t + 1) {
//        while (v.size() % t != 1) {
//            v.push_back(0);
//        }
//
//        vector<int> w2(w.begin(), w.begin() + w.size() - w.size() % t);
//        vector<int> result_temp = _russians(v, w2, s, t, cache);
//        vector<int> b(v.size());
//
//        for (int i = 0; i < b.size(); i++) {
//            b[i] = result_temp[i + 1] - result_temp[i];
//        }
//
//        vector<int> c(w.size() % t, 1);
//        vector<int> d(v.size());
//        vector<int> e(w.size() - w.size() % t);
//        vector<int> w_temp(w.begin() + w.size() - w.size() % t, w.end());
//        vector<int> result(v.size());
//        fill_partial_edit_distance(b, c, v, w_temp, result, e);
//
//        for (int i = 0; i < result.size(); i++) {
//            result[i] += result_temp[0];
//        }
//        return vector(result.begin(), result.begin() + initial_v_size);
//    } else {
//        return edit_distance_rightmost_col(v, w);
//    }
}
