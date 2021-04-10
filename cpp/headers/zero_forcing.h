#ifndef ZERO_FORCING_H
#define ZERO_FORCING_H
#include <vector>
#include <set>
using namespace std;

/* Graph Class */
class graph
{
protected:
    int order;
    vector<pair<int, int>> edges;
    vector< vector<int>> neighbors;
    vector<string> colors;

public:
    int max_degree;
    int min_degree;

    graph() {}

    graph(int order0, const vector<pair<int, int>>& edges0)
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
            vector<int> nbhd;
            for (int i = 0; i < edges.size(); i++) {
                if (edges[i].first == j) {
                    nbhd.push_back(edges[i].second);
                }
                else if (edges[i].second == j) {
                    nbhd.push_back(edges[i].first);
                }
            }

            if (nbhd.size() > max_degree)
                max_degree = nbhd.size();
            if (nbhd.size() < min_degree)
                min_degree = nbhd.size();

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
        vector<int> ret;
        for (int i = 0; i < order; i++)
            ret.push_back(i);

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

    vector<int> adj(int node) const{
        return neighbors[node];
    }

    int get_degree(int node) const {
        return neighbors[node].size();
    }

    void setColor(int node, const string& color) {
        colors[node] = color;
    }

    string getColor(int node) const {
        return colors[node];
    }

    void setAllColor(const string& color) {
        for (int i = 0; i < order; i++)
            colors[i] = color;
    }
	
	void print_nodes()
	{
		cout << "nodes: ";
		for(int i=0; i<order; i++)
		{
			cout << i << " ";
		}
		cout << endl;
	}
	
	void print_edges()
	{
		cout << "edges: ";
		for(int i=0; i<edges.size(); i++)
		{
			cout << "(" << edges[i].first << "," << edges[i].second << ")" << " ";
		}
		cout << endl;
	}
};
/* read graph6 */
graph read_graph6(const string line);
/* Zero Forcing IP */
int zf_ip(graph& g);
/* Zero Forcing Wavefront */
int zf_wave(graph& g);
/* Zero Forcing Heuristic */
set<int,less<int>> heuristic(graph& g);
// returnTriplet Structure
struct returnTriplet {
    int zero_forcing_num;
    int propagation;
    vector<int> zero_forcing_set;
    int throttling_num;
};
/* Zero Forcing GA */
returnTriplet zf_ga(graph& G);

#endif