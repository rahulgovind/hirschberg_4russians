//
// Created by Rahul Govind on 06/12/18.
//

#include "hirscherg.h"
#include "russians.h"

#include<bits/stdc++.h>

using namespace std;

void printResult(map<int, vector<int> > &result) {
    for (auto iter = result.begin(); iter != result.end(); ++iter) {
        cout << iter->first << ": ";
        for (int i = 0; i < iter->second.size(); ++i)
            cout << iter->second[i] << ",";
        cout << endl;
    }
}

void printVector(vector<int> &vec) {
    for (int i = 0; i < vec.size(); ++i)
        cout << vec[i] << ",";
    cout << endl;
}

vector<int> russians_prefix(int i1, int j1, int i2, int j2,
                            string &seq1, string &seq2, int *cache, int t, int s) {
    cout << "Prefix: " << i1 << " " << i2 << " " << j1 << " " << j2 << endl;
    vector<int> v, w;
    for (int i = i1; i <= i2; i++) {
        v.push_back(seq1[i] - '0');
    }
    for (int j = j1; j <= j2; j++) {
        w.push_back(seq2[j] - '0');
    }

    vector<int> result = russians(v, w, s, t, cache);
    for (int i = 0; i < result.size(); i++) {
        result[i] = -result[i];
    }
    return result;
}

vector<int> russians_suffix(int i1, int j1, int i2, int j2,
                            string &seq1, string &seq2, int *cache, int t, int s) {
    cout << "Suffix: " << i1 << " " << i2 << " " << j1 << " " << j2 << endl;
    vector<int> v, w;

    for (int i = i1 + 1; i <= i2; i++) {
        v.push_back(seq1[i] - '0');
    }
    v.push_back(0);
    reverse(v.begin(), v.end());

    for (int j = j1 + 1; j <= j2; j++) {
        w.push_back(seq2[j] - '0');
    }
    w.push_back(0);
    reverse(w.begin(), w.end());

    vector<int> result = russians(v, w, s, t, cache);
    reverse(result.begin(), result.end());
    for (int i = 0; i < result.size(); i++) {
        result[i] = -result[i];
    }
    return result;
}


vector<int> prefix(int i1, int j1, int i2, int j2, string &seq1, string &seq2, map<string, int> &scoring) {
    cout << "prefix: " << i1 << ", " << j1 << ", " << i2 << ", " << j2 << endl;
    // 2-column solution
    int col_length = i2 - i1 + 1;
    vector<int> col1, col2(col_length);
    for (int i = 0; i < col_length; ++i)
        col1.push_back(-i);
    for (int j = j1 + 1; j <= j2; j++) {
        for (int i = 0; i < col_length; ++i) {
            int max = col1[i] + scoring["indel"];
            if (i > 0) { // j > j1 always
                max = max > (col2[i - 1] + scoring["indel"]) ? max : (col2[i - 1] + scoring["indel"]);
                int character_compare_score = 0;
                if (seq1[i1 + i] == seq2[j])
                    character_compare_score = scoring["match"];
                else
                    character_compare_score = scoring["mismatch"];
                max = max > (col1[i - 1] + character_compare_score) ? max : (col1[i - 1] + character_compare_score);
            }
            col2[i] = max;
        }
        if (j < j2)
            col1 = col2;
    }
    printVector(col2);
    return col2;
}

vector<int>  suffix(int i1, int j1, int i2, int j2, string &seq1, string &seq2, map<string, int> &scoring) {
    cout << "suffix: " << i1 << ", " << j1 << ", " << i2 << ", " << j2 << endl;
    // 2-column solution
    int col_length = i2 - i1 + 1;
    vector<int> col1(col_length), col2(col_length);
    for (int i = col_length - 1; i >= 0; --i)
        col1[col_length - 1 - i] = -i;
    for (int j = j2 - 1; j >= j1; --j) {
        for (int i = col_length - 1; i >= 0; --i) {
            int max = col1[i] + scoring["indel"];
            if (i < col_length - 1) {
                max = max > (col2[i + 1] + scoring["indel"]) ? max : (col2[i + 1] + scoring["indel"]);
                int character_compare_score = 0;
                if (seq1[i1 + i + 1] == seq2[j + 1])
                    character_compare_score = scoring["match"];
                else
                    character_compare_score = scoring["mismatch"];
                max = max > (col1[i + 1] + character_compare_score) ? max : (col1[i + 1] + character_compare_score);
            }
            col2[i] = max;
        }
        col1 = col2;
    }
    printVector(col2);
    return col2;
}

void hirschberg(int i1, int j1, int i2, int j2, string &seq1, string &seq2, map<string, int> &scoring,
                map<int, vector<int> > &result, int *cache, int s, int t) {
    cout << "hirsch: " << i1 << ", " << j1 << ", " << i2 << ", " << j2 << endl;
    if (j1 >= j2 - 1)
        return;

//    vector<int> pre = russians_prefix(i1, j1, i2, j1 + ((j2 - j1) / 2), seq1, seq2, cache, t, s);
    vector<int> pre = prefix(i1, j1, i2, j1 + ((j2 - j1) / 2), seq1, seq2, scoring);
//    vector<int> suf = russians_suffix(i1, j1 + ((j2 - j1) / 2), i2, j2, seq1, seq2, cache, t, s);
    vector<int> suf = suffix(i1, j1 + ((j2 - j1) / 2), i2, j2, seq1, seq2, scoring);
//    cout << "SUF1: "; print_vector(suf);
//    cout << "SUF2: "; print_vector(suf2);
    int max = pre[0] + suf[0];
    vector<int> i_star_list;
    i_star_list.push_back(i1);
    for (int i = 1; i < pre.size(); ++i) {
        if (max < pre[i] + suf[i]) {
            max = pre[i] + suf[i];
            i_star_list.clear();
            i_star_list.push_back(i + i1);
        } else if (max == pre[i] + suf[i]) {
            i_star_list.push_back(i + i1);
        }
    }
    int i_star = i_star_list[0];
    result[(j1 + ((j2 - j1) / 2))] = i_star_list;
    printResult(result);
    hirschberg(i1, j1, i_star, (j1 + ((j2 - j1) / 2)), seq1, seq2, scoring, result, cache, s, t);
    hirschberg(i_star, (j1 + ((j2 - j1) / 2)), i2, j2, seq1, seq2, scoring, result, cache, s, t);
}

void addFirstLastColumnToResult(map<int, vector<int> > &result, string &seq1, string &seq2, map<string, int> &scoring) {
    cout << "in addFirstLastColTOResult:" << endl;
    vector<int> i_star_first_col, i_star_last_col;
    vector<int> first_col_suffix = suffix(0, 0, seq1.size() - 1, seq2.size() - 1, seq1, seq2, scoring);
    vector<int> last_col_prefix = prefix(0, 0, seq1.size() - 1, seq2.size() - 1, seq1, seq2, scoring);
    int first_col_best_score = first_col_suffix[0], last_col_best_score = last_col_prefix.back();
    cout << first_col_best_score << " " << last_col_best_score << endl;
    for (int i = 0; i < seq1.size(); ++i) {
        int first_col_score = first_col_suffix[i] - i;
        int last_col_score = last_col_prefix[seq1.size() - i - 1] - i;
        if (i == 0 || first_col_score == first_col_best_score) {
            i_star_first_col.push_back(i);
        }
        if (i == 0 || last_col_score == last_col_best_score) {
            i_star_last_col.push_back(seq1.size() - i - 1);
        }
    }
    reverse(i_star_last_col.begin(), i_star_last_col.end());
    result[0] = i_star_first_col;
    result[seq2.size() - 1] = i_star_last_col;
    cout << "first col best indices:";
    printVector(i_star_first_col);
    cout << "last col best indices:";
    printVector(i_star_last_col);
}

void printAlignment(map<int, vector<int> > &result, string &seq1, string &seq2, map<string, int> &scoring) {
    string seq1_align, seq2_align;
    addFirstLastColumnToResult(result, seq1, seq2, scoring);
    int next_col_first_cell = seq1.size();
    map<int, vector<int> >::reverse_iterator riter = result.rbegin();
    for (; riter != result.rend(); ++riter) {
        int col_no = riter->first;
        bool col_started = false;
        for (int i = riter->second.size() - 1; i >= 0; --i) {
            if (!col_started) {
                if (riter->second[i] > next_col_first_cell)
                    continue;
                if (riter->second[i] == next_col_first_cell) {
                    seq1_align.push_back('-');
                    seq2_align.push_back(seq2[col_no + 1]);
                } else if (riter->second[i] == next_col_first_cell - 1) {
                    if (riter != result.rbegin()) {
                        seq1_align.push_back(seq1[next_col_first_cell]);
                        seq2_align.push_back(seq2[col_no + 1]);
                    }
                }
                col_started = true;
            } else {
                seq1_align.push_back(seq1[riter->second[i + 1]]);
                seq2_align.push_back('-');
            }
        }
        next_col_first_cell = riter->second[0];
        cout << "next_col_first_cell: " << next_col_first_cell << endl;
    }
    reverse(seq1_align.begin(), seq1_align.end());
    reverse(seq2_align.begin(), seq2_align.end());
    cout << seq1_align << endl;
    cout << seq2_align << endl;
}

//int main() {
//    string sequence1, sequence2;
//    cin >> sequence1 >> sequence2;
//    map<string, int> scoring;
//    scoring["indel"] = -1;
//    scoring["mismatch"] = -1;
//    scoring["match"] = 0;
//    sequence1 = "*" + sequence1;
//    sequence2 = "*" + sequence2;
//    if (sequence2.size() < sequence1.size()) // keeping sequence1 as the shorter sequence
//        swap(sequence1, sequence2);
//    cout << "seq1: " << sequence1 << endl;
//    cout << "seq2: " << sequence2 << endl;
//    map<int, vector<int> > result;
//    hirschberg(0, 0, sequence1.size() - 1, sequence2.size() - 1, sequence1, sequence2, scoring, result);
//    for (map<int, vector<int> >::iterator iter = result.begin(); iter != result.end(); ++iter) {
//        cout << "for column " << iter->first << endl;
//        for (int i = 0; i < iter->second.size(); ++i) {
//            cout << iter->second[i] << ", ";
//        }
//        cout << endl;
//    }
//    printAlignment(result, sequence1, sequence2, scoring);
//}