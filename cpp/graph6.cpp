#include <graph6.hpp>
using namespace std;
using namespace wf;
using namespace ip;
/* read graph6 for wavefront */
vector<wf::node> read_graph6_wf(const string &line)
{
	int m = line.size();		// length of line
	int data[m];				// data for graph
	// shift line characters by 63
	for(int i=0; i<m; i++)
	{
		data[i] = line[i] - 63;	
	}
	// order and edge data
	int n;
	int* edge;
	if(data[0]<=62)
	{
		n = data[0];
		edge = &data[1];
		m -= 1;
	}
	else if(data[1]<=62)
	{
		n = (data[1]<<12) + (data[2]<<6) + data[3];
		edge = &data[4];
		m -= 4;
	}
	else
	{
		n = (data[2]<<30) + (data[3]<<24) + (data[4]<<18) + (data[5]<<12) + (data[6]<<6) + data[7];
		edge = &data[8];
		m -= 8;
	}
	// bits
	bool bits[m*6];
	for(int i=0;i<m;i++)
	{
		for(int j=5; j>=0; j--)
		{
			bits[i*6+(5-j)] = (edge[i]>>j) & 1;
		}
	}
	// build vector of nodes
	wf::node blank_node;
	vector<wf::node> nodes(n,blank_node);
	int k = 0;
	for(int j=1; j<n; j++)
	{
		for(int i=0; i<j; i++)
		{
			if(bits[k])
			{
				nodes[i].add_to_adj_list(j);
				nodes[j].add_to_adj_list(i);
			}
			k++;
		}
	}
	// return
	return nodes;
}
/* read graph6 for fort_cover_ip */
ip::graph read_graph6_ip(const string &line)
{
	int m = line.size();		// length of line
	int data[m];				// data for graph
	// shift line characters by 63
	for(int i=0; i<m; i++)
	{
		data[i] = line[i] - 63;	
	}
	// order and edge data
	int n;
	int* edge;
	if(data[0]<=62)
	{
		n = data[0];
		edge = &data[1];
		m -= 1;
	}
	else if(data[1]<=62)
	{
		n = (data[1]<<12) + (data[2]<<6) + data[3];
		edge = &data[4];
		m -= 4;
	}
	else
	{
		n = (data[2]<<30) + (data[3]<<24) + (data[4]<<18) + (data[5]<<12) + (data[6]<<6) + data[7];
		edge = &data[8];
		m -= 8;
	}
	// bits
	bool bits[m*6];
	for(int i=0;i<m;i++)
	{
		for(int j=5; j>=0; j--)
		{
			bits[i*6+(5-j)] = (edge[i]>>j) & 1;
		}
	}
	// build graph
	fcip::graph our_graph;
	fcip::node blank_node;
	vector<fcip::node> nodes(n,blank_node);
	our_graph.nodes = nodes;
	int k = 0;
	for(int j=1; j<n; j++)
	{
		for(int i=0; i<j; i++)
		{
			if(bits[k])
			{
				our_graph.edges.add_edge(i,j);
				our_graph.edges.add_edge(j,i);
				our_graph.nodes[i].add_to_adj_list(j);
				our_graph.nodes[j].add_to_adj_list(i);
			}
			k++;
		}
	}
	// return graph
	return our_graph;
}
/* read graph6 for genetic_algorithm*/
Graph read_graph6_ga(const string &line)
{
	int m = line.size();		// length of line
	int data[m];				// data for graph
	// shift line characters by 63
	for(int i=0; i<m; i++)
	{
		data[i] = line[i] - 63;	
	}
	// size and edge data
	int n;
	int* edge;
	if(data[0]<=62)
	{
		n = data[0];
		edge = &data[1];
		m -= 1;
	}
	else if(data[1]<=62)
	{
		n = (data[1]<<12) + (data[2]<<6) + data[3];
		edge = &data[4];
		m -= 4;
	}
	else
	{
		n = (data[2]<<30) + (data[3]<<24) + (data[4]<<18) + (data[5]<<12) + (data[6]<<6) + data[7];
		edge = &data[8];
		m -= 8;
	}
	// bits
	bool bits[m*6];
	for(int i=0;i<m;i++)
	{
		for(int j=5; j>=0; j--)
		{
			bits[i*6+(5-j)] = (edge[i]>>j) & 1;
		}
	}
	// build graph edges
	vector<pair<int,int>> edges;
	int k = 0;
	for(int j=1; j<n; j++)
	{
		for(int i=0; i<j; i++)
		{
			if(bits[k])
			{
				edges.push_back(make_pair(i,j));
			}
			k++;
		}
	}
	
	// return graph
	return Graph(n,edges);
}