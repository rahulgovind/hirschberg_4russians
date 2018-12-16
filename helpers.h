//
// Created by Rahul Govind on 15/12/18.
//

#ifndef PROJECT_HELPERS_H
#define PROJECT_HELPERS_H

#include <bits/stdc++.h>

using namespace std;

class Encoder {
private:
    map<char, char> e;
    map<char, char> d;
public:
    Encoder(string s) {
        char present[256];
        for (int i = 0; i < 256; i++) { present[i] = false; }
        for (int i = 0; i < s.length(); i++) { present[s[i]] = true; }

        int count = 0;
        for (int i = 0; i < 256; i++) {
            if (present[i]) {
                e[(char) i] = (char) ('0' + count);
                d[(char) ('0' + count)] = (char) i;
                count++;
            }
        }
    }

    size_t charset_size() {
        return e.size();
    }

    string encode(string s) {
        vector<char> v;
        for (int i = 0; i < s.length(); i++) {
            v.push_back(e[s[i]]);
        }
        return string(v.begin(), v.end());
    }

    string decode(string s) {
        vector<char> v;
        for (int i = 0; i < s.length(); i++) {
            v.push_back(d[s[i]]);
        }
        return string(v.begin(), v.end());
    }
};

#endif //PROJECT_HELPERS_H
