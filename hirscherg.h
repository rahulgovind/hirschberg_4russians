//
// Created by Rahul Govind on 06/12/18.
//

#ifndef PROJECT_HIRSCHERG_H
#define PROJECT_HIRSCHERG_H

#include <bits/stdc++.h>

using namespace std;


typedef pair<string, string> pstring;


pstring hirschberg_standard(string seq1, string seq2);

pstring hirschberg_standard(string seq1, string seq2, map<string, int> &scoring);

pstring hirschberg_russians(string seq1, string seq2, int s, int t);

void printAlignment(map<int, vector<int> > &result, string &seq1, string &seq2, map<string, int> &scoring);

#endif //PROJECT_HIRSCHERG_H
