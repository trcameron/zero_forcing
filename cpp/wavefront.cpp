#include <wavefront.hpp>
using namespace std;
using namespace wf;
int wf::zero_forcing(std::vector<wf::node> nodes)
{
	/******Wavefront (dynamic programming) forcing***************************************/
	int num_nodes = nodes.size();
    std::unordered_map<size_t, wf::node_subset > closures;

    int max_degree = -1;
    for(int i=0; i<num_nodes; ++i)
	{
		if(nodes[i].get_degree() > max_degree)
		{
			max_degree = nodes[i].get_degree();
		}
	}
	if(max_degree==0)
	{
		return num_nodes;
	}
	
	size_t best_set =0;
	wf::node_subset a_subset;
    a_subset.set_weight();
    a_subset.update_closure(nodes);
	closures.insert(std::make_pair(a_subset.create_hash(),a_subset));

    int stop_criteria = 0;
    std::vector<wf::node_subset>::iterator best_subset_iterator;
    int best_size = num_nodes+5;
	
	for(int R= 1; R<num_nodes; ++R)
	{
		std::unordered_map<size_t, wf::node_subset> temp_map(closures);
		for(std::unordered_map<size_t, wf::node_subset>::iterator it=closures.begin(); it!=closures.end(); ++it)
		{
			for(int j=0; j<num_nodes; ++j)
			{
				int unfilled_neighbor = -1;
				const int k = (it->second).num_unfilled_neighbors(nodes[j], unfilled_neighbor);
				//Case for j in closure*************************************************
				if((it->second).is_in_closure(j) && k-1 <= R - (it->second).get_weight())
				{
					wf::node_subset new_subset;
					for(std::set<int>::iterator jt = (it->second).get_subset_begin(); jt != (it->second).get_subset_end(); ++jt)
					{
						new_subset.add_to_subset(*jt);
					}
					for(int l=0; l<nodes[j].get_degree(); ++l)
					{
						const int neighbor = nodes[j].get_adj_list_member(l);
						if(neighbor != unfilled_neighbor && !(it->second).is_in_closure(neighbor))
						{
							new_subset.add_to_subset(neighbor);
						}
					}
					new_subset.set_weight();
	      		  	new_subset.update_closure(nodes);
					
					const size_t new_hash = new_subset.create_hash();
					if(closures.count(new_hash) < 1)
					{
						temp_map.insert(std::make_pair(new_hash, new_subset));
					}
				}
				//Case for j not in closure*************************************************
	    		if(!(it->second).is_in_closure(j) && k <= R - (it->second).get_weight())
				{
					wf::node_subset new_subset;
					for(std::set<int>::iterator jt = (it->second).get_subset_begin(); jt != (it->second).get_subset_end(); ++jt)
					{
						new_subset.add_to_subset(*jt);
					}
					new_subset.add_to_subset(j);
	      		  	for(int l=0; l<nodes[j].get_degree(); ++l)
					{
						const int neighbor = nodes[j].get_adj_list_member(l);
						if(neighbor != unfilled_neighbor && !(it->second).is_in_closure(neighbor))
						{
							new_subset.add_to_subset(neighbor);
						}
					}
					new_subset.set_weight();
					new_subset.update_closure(nodes);
	      		  	const size_t new_hash = new_subset.create_hash();

	      		  	if(closures.count(new_hash) < 1)
					{
						temp_map.insert(std::make_pair(new_hash, new_subset));
					}
				}
			}
			if(it->second.get_weight() < R - max_degree)
			{
				temp_map.erase(it->first);
			}
		}
		closures.swap(temp_map);
		//Check if one of closures is entire graph
		for(std::unordered_map<size_t, wf::node_subset>::iterator it = closures.begin(); it != closures.end(); ++it)
		{
			if((it->second).get_closure_size() == num_nodes && (it->second).get_weight() < best_size)
			{
				stop_criteria = 1;
				best_size = (it->second).get_weight();
				best_set = it->first;
			}
		}
		if(stop_criteria == 1)
		{
			break;
		}
	}
	// return zero forcing number
	return closures[best_set].get_weight();
}