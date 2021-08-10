#include <edg.hpp>
#include <iostream>
#include <fstream>
using namespace std;
using namespace wf;
using namespace ip;
using namespace ga;
/* read edg for wavefront */
std::vector<wf::node> read_edg_wf(const std::string &file_name)
{
	int m, n; 				// number of edges and nodes
	fstream file; 			// edg file
	file.open(file_name);
    if(file.is_open())
	{
		file >> n;
		file >> m;
	}
    else
	{
		cerr << "ERROR: could not open file: " << file_name << endl;
    }
	wf::node blank_node;
	vector<wf::node> nodes(n,blank_node);
	int i,j,k;
	for(k=0; k<m; k++)
	{
		file >> i;
		file >> j;
		nodes[i].add_to_adj_list(j);
		nodes[j].add_to_adj_list(i);
	}
	file.close();
	return nodes;
}
/* read edg for fort_cover_ip */
ip::graph read_edg_ip(const std::string &file_name)
{
	int m, n; 				// number of edges and nodes
	fstream file; 			// edg file
	file.open(file_name);
    if(file.is_open())
	{
		file >> n;
		file >> m;
	}
    else
	{
		cerr << "ERROR: could not open file: " << file_name << endl;
    }
	ip::graph our_graph;
	ip::node blank_node;
	vector<ip::node> nodes(n,blank_node);
	our_graph.nodes = nodes;
	int i,j,k;
	for(k=0; k<m; k++)
	{
		file >> i;
		file >> j;
		nodes[i].add_to_adj_list(j);
		nodes[i].add_to_adj_list(i);
		our_graph.edges.add_edge(i,j);
		our_graph.edges.add_edge(j,i);
		our_graph.nodes[i].add_to_adj_list(j);
		our_graph.nodes[j].add_to_adj_list(i);
	}
	file.close();
	return our_graph;
}
/* read edg for genetic_algorithm */
ga::Graph read_edg_ga(const std::string &file_name)
{
	int m, n; 				// number of edges and nodes
	fstream file; 			// edg file
	file.open(file_name);
    if(file.is_open())
	{
		file >> n;
		file >> m;
	}
    else
	{
		cerr << "ERROR: could not open file: " << file_name << endl;
    }
	vector<pair<int,int>> edges;
	int i,j,k;
	for(k=0; k<m; k++)
	{
		file >> i;
		file >> j;
		edges.push_back(make_pair(i,j));
	}
	file.close();
	return Graph(n,edges);
}