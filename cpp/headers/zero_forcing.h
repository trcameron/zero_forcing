#ifndef ZERO_FORCING_H
#define ZERO_FORCING_H
#include <vector>
using namespace std;

/* Graph Class */
class graph
{
	protected:
		int order;
		vector<pair<int,int>> edges;
		
	public:
		graph(){}
		graph(int order0)
		{
			order = order0;
		}
		graph(int order0,vector<pair<int,int>> edges0)
		{
			order = order0;
			for(int i=0; i<edges0.size(); i++)
			{
				edges.push_back(edges0[i]);
			}
		}
		void add_edge(pair<int,int> e)
		{
			edges.push_back(e);
		}
		int get_order()
		{
			return order;
		}
		vector<pair<int,int>> get_edges()
		{
			return edges;
		}
		vector<int> get_neighbors(int node)
		{
			vector<int> nbhd;
			for(int i=0;i<edges.size();i++)
			{
				if(edges[i].first==node)
				{
					nbhd.push_back(edges[i].second);
				}
				else if(edges[i].second==node)
				{
					nbhd.push_back(edges[i].first);
				}
			}
			return nbhd;
		}
		int get_degree(int node)
		{
			int deg = 0;
			for(int i=0; i<edges.size(); i++)
			{
				if(edges[i].first==node || edges[i].second==node)
				{
					deg += 1;
				}
			}
			return deg;
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

graph read_graph6(const string line);
int zf_ip(graph g);

#endif