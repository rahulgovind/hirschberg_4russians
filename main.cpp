#include <iostream>
#include <bits/stdc++.h>
#include <cstdlib>
#include <random>
#include "russians.h"
#include "hirscherg.h"
#include "helpers.h"
#include <chrono>
#include "cxxopts.hpp"

using namespace std;
using namespace std::chrono;

typedef long long ll;

typedef vector<int> vint;
typedef pair<vint, vint> vpair;

int edit_distance2(const vector<int> &v, const vector<int> &w) {
    // One column at a time
    ulong m = v.size() + 1, n = w.size() + 1;
    vector<int> prev(m), curr(m);

    for (int i = 0; i < m; i++) {
        prev[i] = i;
    }

    for (int j = 1; j < n; j++) {
        curr[0] = j;
        for (int i = 1; i < m; i++) {
            if (v[i - 1] == w[j - 1]) {
                curr[i] = prev[i - 1];
            } else {
                curr[i] = min(curr[i - 1], min(prev[i], prev[i - 1])) + 1;
            }
        }
        prev = curr;
    }

    return curr.back();
}

int edit_distance(const vector<int> &v, const vector<int> &w) {
    unsigned long m = v.size() + 1, n = w.size() + 1;
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
            if (v[i - 1] == w[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = min(dp[i - 1][j], min(dp[i][j - 1], dp[i - 1][j - 1])) + 1;
            }
        }
    }

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
    uniform_int_distribution<mt19937::result_type> dist(0, s - 1);

    vector<int> result;
    for (int i = 0; i < n; i++) {
        result.push_back(dist(rng));
    }
    return result;
}

string generate_random_string(int n, int s) {
    mt19937 rng;
    rng.seed(std::random_device()());
    uniform_int_distribution<mt19937::result_type> dist(0, s - 1);

    vector<char> result;
    for (int i = 0; i < n; i++) {
        result.push_back('0' + (char) dist(rng));
    }
    return string(result.begin(), result.end());
}

int char_count(string s) {
    int r = 0;
    for (int i = 0; i < s.length(); i++) {
        if (s[i] != '-')
            r++;
    }
    return r;
}

void print_alignment(pstring alignment) {
    cout << alignment.first << endl;
    cout << alignment.second << endl;
}

// Runs hirschberg_russians against some user-defined strings
void test1() {
    string sequence1 = "001011011001001";
    for (int i = 0; i < 3; i++) {
        sequence1 += sequence1;
    }
    string sequence2 = "001011110101001";
    for (int i = 0; i < 6; i++) {
        sequence2 += sequence2;
    }
    map<int, vector<int> > result;

    pstring alignment = hirschberg_russians(sequence1, sequence2, 4, 2);
    print_alignment(alignment);
    cerr << "Edit distance naive: " << edit_distance(string_to_vector(sequence1), string_to_vector(sequence2)) << endl;
}

// Runs hirschberg against two randomly generated strings and times them
void test2() {
    string seq1 = generate_random_string(10000, 2);
    string seq2 = generate_random_string(10000, 2);
    cerr << "Sequence 1: " << seq1 << endl;
    cerr << "Sequence 2: " << seq2 << endl;


    auto t1 = high_resolution_clock::now();
    pstring alignment = hirschberg_russians(seq1, seq2, 2, 3);
    auto t2 = high_resolution_clock::now();
    printf("Hirschberg (Four Russians) took %lld milliseconds\n", duration_cast<milliseconds>(t2 - t1).count());

    t1 = high_resolution_clock::now();
    alignment = hirschberg_standard(seq1, seq2);
    t2 = high_resolution_clock::now();
    printf("Hirschberg (Standard) took %lld milliseconds\n", duration_cast<milliseconds>(t2 - t1).count());

    t1 = high_resolution_clock::now();
    edit_distance(string_to_vector(seq1), string_to_vector(seq2));
    t2 = high_resolution_clock::now();
    printf("Naive Edit Distance took %lld milliseconds\n", duration_cast<milliseconds>(t2 - t1).count());

    cerr << "Edit distance naive2:\t" << edit_distance2(string_to_vector(seq1), string_to_vector(seq2)) << endl;
}

// Reads from input
void test3() {
    string s1, s2;
    cin >> s1;
    cin >> s2;

    Encoder e(s1 + s2);
    int s = (int) e.charset_size();
    string seq1 = e.encode(s1);
    string seq2 = e.encode(s2);

    cerr << "Sequence 1: " << seq1 << endl;
    cerr << "Sequence 2: " << seq2 << endl;

    auto t1 = high_resolution_clock::now();
    pstring alignment = hirschberg_russians(seq1, seq2, s, 3);
    auto t2 = high_resolution_clock::now();
    printf("Hirschberg (Four Russians) took %lld milliseconds\n", duration_cast<milliseconds>(t2 - t1).count());

    t1 = high_resolution_clock::now();
    alignment = hirschberg_standard(seq1, seq2);
    t2 = high_resolution_clock::now();
    printf("Hirschberg (Standard) took %lld milliseconds\n", duration_cast<milliseconds>(t2 - t1).count());

    t1 = high_resolution_clock::now();
    edit_distance(string_to_vector(seq1), string_to_vector(seq2));
    t2 = high_resolution_clock::now();
    printf("Naive Edit Distance took %lld milliseconds\n", duration_cast<milliseconds>(t2 - t1).count());

    cerr << "Edit distance naive2:\t" << edit_distance2(string_to_vector(seq1), string_to_vector(seq2)) << endl;

}

void cli(int argc, char **argv) {
    cxxopts::Options options("Hirschberg-4Russians",
                             "Optimal edit distance algorithm with hirschberg and four russians technique");
    options.add_options()("m,method", "Method to use", cxxopts::value<std::string>());
    options.add_options()("f,file", "Input file", cxxopts::value<std::string>());

    auto result = options.parse(argc, argv);
    string method = result["method"].as<std::string>();

    string s1, s2;
    
    if (result.count("file") > 0) {
        string fname = result["file"].as<std::string>();
        cerr << "Reading from file: " << fname << endl;

        ifstream f(fname);
        f >> s1;
        f >> s2;
    } else {
        cin >> s1;
        cin >> s2;
    }

    Encoder e(s1 + s2);
    int s = (int) e.charset_size();
    string seq1 = e.encode(s1);
    string seq2 = e.encode(s2);

    cerr << "Sequence 1: " << seq1 << endl;
    cerr << "Sequence 2: " << seq2 << endl;

    if (method == "russians") {
        auto t1 = high_resolution_clock::now();
        pstring alignment = hirschberg_russians(seq1, seq2, s, 3);
        auto t2 = high_resolution_clock::now();
        printf("Hirschberg (Four Russians) took %lld milliseconds\n", duration_cast<milliseconds>(t2 - t1).count());
    } else if (method == "hirschberg") {
        auto t1 = high_resolution_clock::now();
        auto alignment = hirschberg_standard(seq1, seq2);
        auto t2 = high_resolution_clock::now();
        printf("Hirschberg (Standard) took %lld milliseconds\n", duration_cast<milliseconds>(t2 - t1).count());
    } else {
        auto t1 = high_resolution_clock::now();
        edit_distance(string_to_vector(seq1), string_to_vector(seq2));
        auto t2 = high_resolution_clock::now();
        printf("Naive Edit Distance took %lld milliseconds\n", duration_cast<milliseconds>(t2 - t1).count());

        cerr << "Edit distance naive2:\t" << edit_distance2(string_to_vector(seq1), string_to_vector(seq2)) << endl;
    }
}


int main(int argc, char **argv) {
    cli(argc, argv);
}