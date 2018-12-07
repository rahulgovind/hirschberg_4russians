#include <iostream>
#include <bits/stdc++.h>
#include <cstdlib>
#include <random>
#include "russians.h"
#include "hirscherg.h"

using namespace std;

typedef long long ll;

typedef vector<int> vint;
typedef pair<vint, vint> vpair;



int edit_distance(const vector<int> &v, const vector<int> &w) {
    unsigned long m = v.size(), n = w.size();
    vector<vector<int> > dp(m, vector<int>(n));

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

//    print_matrix(dp);
    return dp[m - 1][n - 1];
}




vector<int> string_to_vector(const string &s) {
    vector<int> result(s.length());
    for (int i = 0; i < s.length(); i++) {
        result[i] = s[i];
    }
    return result;
}

vector<int> generate_random_array(int n, int s) {
    mt19937 rng;
    rng.seed(std::random_device()());
    uniform_int_distribution<mt19937::result_type> dist(0, s-1);

    vector<int> result;
    for (int i=0; i < n; i++) {
        result.push_back(dist(rng));
    }
    return result;
}

int main() {
    int t = 2;
    int s = 4;
    int *cache = cache_values(t, s);

//    string sequence1, sequence2;
//    cin >> sequence1 >> sequence2;
    string sequence1 = "001011011001001";
    for (int i=0; i < 3; i++) {
        sequence1 += sequence1;
    }
    string sequence2 = "001011110101001";

    for (int i=0; i < 3; i++) {
        sequence2 += sequence2;
    }
    map<string, int> scoring;
    scoring["indel"] = -1;
    scoring["mismatch"] = -1;
    scoring["match"] = 0;
    sequence1 = sequence1;
    sequence2 = sequence2;
    if (sequence2.size() < sequence1.size()) // keeping sequence1 as the shorter sequence
        swap(sequence1, sequence2);
    cout << "seq1: " << sequence1 << endl;
    cout << "seq2: " << sequence2 << endl;
    map<int, vector<int> > result;
    hirschberg(0, 0, sequence1.size() - 1, sequence2.size() - 1, sequence1, sequence2, scoring, result, cache, s, t);
    for (map<int, vector<int> >::iterator iter = result.begin(); iter != result.end(); ++iter) {
        cout << "for column " << iter->first << endl;
        for (int i = 0; i < iter->second.size(); ++i) {
            cout << iter->second[i] << ", ";
        }
        cout << endl;
    }
    printAlignment(result, sequence1, sequence2, scoring);
}