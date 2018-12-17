//
// Created by Rahul Govind on 05/12/18.
//

#ifndef PROJECT_RUSSIANS_H
#define PROJECT_RUSSIANS_H

#include <vector>
#include <map>

typedef unsigned long long ulong;

void print_vector(std::vector<int> v);

std::vector<int> russians(std::vector<char> v, std::vector<char> w, int s, int t, const char *cache);

char *cache_values(int t, int s);

char *calculate_or_load_cache(int t, int s);

#endif //PROJECT_RUSSIANS_H
