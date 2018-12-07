//
// Created by Rahul Govind on 06/12/18.
//

#ifndef PROJECT_HIRSCHERG_H
#define PROJECT_HIRSCHERG_H

#include <bits/stdc++.h>

using namespace std;

void hirschberg(int i1, int j1, int i2, int j2, string &seq1, string &seq2, map<string, int> &scoring,
                map<int, vector<int> > &result, int *cache, int s, int t);

void printAlignment(map<int, vector<int> > &result, string &seq1, string &seq2, map<string, int> &scoring);

#endif //PROJECT_HIRSCHERG_H
