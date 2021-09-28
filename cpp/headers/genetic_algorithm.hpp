#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <stdlib.h>
#include <utility>
#include <list>
#include <fstream>
#include <bitset>
#include <unordered_map>
using namespace std;
namespace ga
{
/* constant expression */
constexpr auto M = 300;
/* Bitset structure */
struct Bitset {
    bitset<M> bit;
    int size;

    Bitset(){}

    Bitset(int order) { size = order; }

    Bitset(vector<int>& v, int order) {
        size = order;
        for (int i = 0; i < v.size(); i++)
            bit.set(v[i]);
    }

    vector<int> convert() {
        vector<int> output;
        for (int i = 0; i < size; i++)
            if (bit[i])
                output.push_back(i);
        return output;
    }

    bool operator==(const Bitset& rhs) const {
        return bit == rhs.bit && size == rhs.size;
    }

    bool operator!=(const Bitset& rhs) const {
        return !(*this == rhs);
    }

    friend ostream& operator<<(ostream& os, const Bitset& b) {
        cout << "[";
        for (int i = 0; i < b.size; i++) {
            if (!b.bit[i])
                continue;
            cout << i << ", ";
        }
        cout << "]";
        return os;
    }
};
/* Bitset operator */
Bitset operator| (const Bitset& lhs, const Bitset& rhs);
/* sample functions */
void sample(vector<int>& vertices, Bitset& newVector, int k);
void sample(vector<int>& vertices, vector<int>& newVector, int k);
/* Graph class */
class Graph
{
protected:
    int order;
    vector<pair<int, int>> edges;
    vector<bool> colors;

public:
    vector< Bitset> neighbors;
    int max_degree;
    int min_degree;

    Graph() {}

    Graph(int order0, const vector<pair<int, int>>& edges0)
    {
        order = order0;
        for (int i = 0; i < edges0.size(); i++)
        {
            edges.push_back(edges0[i]);
        }

        neighbors = {};
        max_degree = 0;
        min_degree = 9999999;
        for (int j = 0; j < order; j++) {
            colors.push_back("white");
            Bitset nbhd(order);
            for (int i = 0; i < edges.size(); i++) {
                if (edges[i].first == j) {
                    nbhd.bit.set(edges[i].second);
                }
                else if (edges[i].second == j) {
                    nbhd.bit.set(edges[i].first);
                }
            }

            if (nbhd.bit.count() > max_degree)
                max_degree = nbhd.bit.count();
            if (nbhd.bit.count() < min_degree)
                min_degree = nbhd.bit.count();

            neighbors.push_back(nbhd);
        }
    }
    
    int get_order() const
    {
        return order;
    }

    vector<pair<int, int>> get_edges() const
    {
        return edges;
    }

    vector<int> vertices() const {
        vector<int> ret(order, 0);
        for (int i = 0; i < order; i++)
            ret[i] = i;

        return ret;
    }

    vector<vector<int>> get_adj() const
    {
        vector<vector<int>> adj(order, vector<int>(order, 0));
        for (int k = 0; k < edges.size(); k++)
        {
            adj[edges[k].first][edges[k].second] = 1;
            adj[edges[k].second][edges[k].first] = 1;
        }
        return adj;
    }

    Bitset adj(int node) const{
        return neighbors[node];
    }

    int get_degree(int node) const {
        return neighbors[node].bit.count();
    }

    void setColor(int node, const bool& color) {
        colors[node] = color;
    }

    bool getColor(int node) const {
        return colors[node];
    }

    void setAllColor(const bool& color) {
        for (int i = 0; i < order; i++)
            colors[i] = color;
    }
};
/* white_path */
bool white_path(Graph& G, int s, int d);
/* std_rule */
bool std_rule(Graph& G, int node0, const Bitset& b);
/* psd_rule */
bool psd_rule(Graph& G, int node0, const Bitset& b);
/* skew_rule */
bool skew_rule(Graph& G, int node0, const Bitset& b);
/* returnPair structure */
struct returnPair {
    Bitset vertices;
    int propagation;
};
/* forcing_process */
returnPair forcing_process(Graph& G, const Bitset& b, bool (*rule)(Graph&, int, const Bitset&) = &std_rule, int t = 1);
/* heuristic */
Bitset heuristic(Graph&G, bool (*rule)(Graph&, int, const Bitset&) = &std_rule);
/* returnQuad structure */
struct returnQuad {
    int zero_forcing_num;
    int propagation;
    vector<int> zero_forcing_set;
    int throttling_num;
};
/* zero_forcing */
returnQuad zero_forcing(Graph& G, bool (*rule)(Graph&, int, const Bitset&) = &std_rule, bool throttling_num = false);
};
#endif
