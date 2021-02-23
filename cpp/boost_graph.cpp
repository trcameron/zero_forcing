#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost;

int main(int argc,char** argv)
{
	// create a typedef for the Graph type
	typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
	
	// Make convenient labels for the vertices
  	enum { A, B, C, D, E, N };
	const int num_vertices = N;
	const char* name = "ABCDE";
	
	// writing out the edges in the graph
	typedef std::pair<int, int> Edge;
	Edge edge_array[] = {Edge(A,B),Edge(A,D),Edge(C,A),Edge(D,C),Edge(C,E),Edge(B,D),Edge(D,E) };
	const int num_edges = sizeof(edge_array)/sizeof(edge_array[0]);
	
	// declare a graph object
	Graph g(edge_array, edge_array + sizeof(edge_array) / sizeof(Edge), num_vertices);
	
	typedef graph_traits<Graph>::vertex_descriptor Vertex;
	// get the property map for vertex indices
	typedef property_map<Graph, vertex_index_t>::type IndexMap;
	IndexMap index = get(vertex_index, g);
	
	std::cout << "vertices(g) = ";
	typedef graph_traits<Graph>::vertex_iterator vertex_iter;
	std::pair<vertex_iter, vertex_iter> vp;
	for (vp = vertices(g); vp.first != vp.second; ++vp.first)
	{
		Vertex v = *vp.first;
		std::cout << index[v] <<  " ";
	}
	std::cout << std::endl;
	
	std::cout << "edges(g) = ";
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
		std::cout << "(" << index[source(*ei, g)] << "," << index[target(*ei, g)] << ") ";
	    std::cout << std::endl;
	
	return 0;
}