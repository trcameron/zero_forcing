#ifndef ZERO_FORCING_H
#define ZERO_FORCING_H
#include <vector>
using namespace std;

vector<vector<int>> read_graph6(const string line);
int zf_ip(vector<vector<int>> adj);

#endif