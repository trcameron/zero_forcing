#ifndef WAVEFRONT_H
#define WAVEFRONT_H
#include<set>
#include<unordered_map>
#include<vector>

using namespace std;

namespace wf
{
class node{

public:
  int get_degree(){
    return adj_list.size();
  };
  void add_to_adj_list(int i){
    adj_list.push_back(i);
  };
  int get_adj_list_member(int i){
    return adj_list[i];
  };

private:
  std::vector<int> adj_list;

};
class node_subset{

public:
	std::vector<bool> closure_contains;
  int weight;
  std::set<int> closure;
  std::set<int> subset;

  void set_weight(){
    weight = subset.size();
  };
  int get_weight(){
    return weight;
  };
  bool is_in_closure(int i){
    return closure_contains[i];
  };

  void add_to_subset(int i){
    subset.insert(i);
  };

  void add_to_closure(int i){
    closure.insert(i);
  };

  int num_unfilled_neighbors(node& v, int& unfilled_neighbor){

    const int v_degree = v.get_degree();
    int result = 0;
    for(int i=0; i<v_degree; ++i){
      if(!closure_contains[v.get_adj_list_member(i)]){
	result++;
	unfilled_neighbor = v.get_adj_list_member(i);
      }
    }
    return result;
  };

  int get_closure_size(){
    return closure.size();
  }

  std::set<int>::iterator get_closure_begin(){
    return closure.begin();
  };
  std::set<int>::iterator get_closure_end(){
    return closure.end();
  };
  std::set<int>::iterator get_subset_begin(){
    return subset.begin();
  };
  std::set<int>::iterator get_subset_end(){
    return subset.end();
  };

  void update_closure(std::vector<node>& nodes){

    const int num_nodes = nodes.size();
    closure_contains = std::vector<bool>(num_nodes, false);
    for(std::set<int>::iterator it = subset.begin(); it != subset.end(); ++it){
      closure.insert(*it);
      closure_contains[*it] = true;
    }

    bool end_flag = false;
    while(!end_flag){
      end_flag = true;
      std::set<int> temp_set(closure);
      std::vector<bool> temp_closure_contains(closure_contains);
      for(std::set<int>::iterator it = closure.begin(); it != closure.end(); ++it){
	//Count the number of forced nodes that are adjacent to the node in the closure
	int unforced_sum = 0;
	int white_neighbor;
	const int node_degree = nodes[*it].get_degree();
	for(int i=0; i<node_degree; ++i){
	  const int neighbor = nodes[*it].get_adj_list_member(i);
	  if(!closure_contains[neighbor]){
	    unforced_sum++;
	    white_neighbor = neighbor;
	  }
	  if(unforced_sum==2){
	    break;
	  }
	}
	if(unforced_sum == 1){

	  temp_set.insert(white_neighbor);

	  temp_closure_contains[white_neighbor] = true;
	  end_flag = false;
	}
      }
      closure.swap(temp_set);
      closure_contains.swap(temp_closure_contains);
    }

  };

  size_t create_hash(){
    std::hash< std::vector<bool> > hasher;
    return hasher(closure_contains);
  };

private:


};
/* zero_forcing */
int zero_forcing(std::vector<wf::node> nodes);
};
#endif