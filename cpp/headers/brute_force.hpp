#ifndef BRUTE_FORCE_H
#define BRUTE_FORCE_H
#include <boost/functional/hash.hpp>
namespace bf
{
/* type definitions */
typedef pair<int,int> edge;
typedef unordered_set<edge,boost::hash<edge>> edge_set;
typedef unordered_set<int> node_set;
typedef unordered_map<int,node_set> nbhd_map;
/* graph class */
class graph
{
protected:
	int order = 0, size = 0;
	node_set nodes;
	edge_set edges;
	nbhd_map neighbors;
public:
	graph(){}
	graph(const int order0,const edge_set& edges0)
	{
		order = order0;
		edges = edges0;
		size = edges.size();
		for(int i=0; i<order; i++)
		{
			nodes.insert(i);
			node_set temp;
			for(edge_set::iterator it=edges.begin(); it!=edges.end(); it++)
			{
				if(it->first==i)
				{
					temp.insert(it->second);
				}
				else if(it->second==i)
				{
					temp.insert(it->first);
				}
			}
			neighbors[i] = temp;
			colored[i] = false;
		}	
	}
	graph(const node_set& nodes0,const edge_set& edges0)
	{
		nodes = nodes0;
		edges = edges0;
		order = nodes.size();
		size = edges.size();
		for(node_set::iterator nit=nodes.begin(); nit!=nodes.end(); nit++)
		{
			node_set temp;
			for(edge_set::iterator eit=edges.begin(); eit!=edges.end(); eit++)
			{
				if(eit->first==*nit)
				{
					temp.insert(eit->second);
				}
				else if(eit->second==*nit)
				{
					temp.insert(eit->first);
				}
			}
			neighbors[*nit] = temp;
			colored[*nit] = false;
		}
	}
	const node_set get_neighbors(const int u) const
	{
		return neighbors.at(u);
	}
	int get_order() const
	{
		return order;
	}
	int get_size() const
	{
		return size;
	}
};
};
#endif