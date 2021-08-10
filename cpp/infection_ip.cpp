#include <iostream>
#include <fstream>
#include <ip.hpp>
#include <infection_ip.hpp>
#define STOPPING_TIME 7200
#define PROBLEM_TYPE 1 //1 for IP, 0 for LP
using namespace std;
using namespace ip;
using namespace ifip;
/* cut_callback */
int __stdcall ifip::cut_callback(GRBmodel *model, void *cbdata, int where, void *usrdata){


  if(where == GRB_CB_MIPNODE){

    double running_time;
    int error;
      error = GRBcbget(cbdata, where, GRB_CB_RUNTIME, &running_time);

      if(error){
        cerr << "ERROR: " << error << " getting runtime." << endl;
      }
      if(running_time > STOPPING_TIME){
	      GRBterminate(model);
      }
  }

  return 0;
};
/* zero_forcing */
int ifip::zero_forcing(ip::graph &our_graph)
{
    GRBenv *env = NULL;
    GRBmodel *zero_forcing = NULL;
    GRBmodel *zero_forcing_linear = NULL;
    GRBmodel *col_gen = NULL;
    int error = 0;
	int num_nodes = our_graph.nodes.size();
    int num_edges = our_graph.edges.get_size();
	int max_timesteps = num_nodes;

    //find min degree
    int min_degree = num_nodes;
    for(int i=0; i<num_nodes; ++i){
      if(our_graph.nodes[i].get_degree() < min_degree){
	min_degree = our_graph.nodes[i].get_degree();
      }
    }

	//Gurobi model********************************************************************************
    //Create the Gurobi environment
    error = GRBloadenv(&env, NULL);
    if(error){
      cout << "ERROR at: load environment " << error << endl;
      return 1;
    }
    //error = GRBsetintparam(GRBgetenv(zero_forcing), "UpdateMode", 1);
    //const int num_variables = num_nodes*3+num_edges;


    /*Create the model*******************************************************/
    const int num_variables = 2*num_nodes + num_edges;
    const int s_start = 0;
    const int x_start = s_start + num_nodes;
    const int y_start = x_start + num_nodes;

    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    for(int i=0; i<x_start; ++ i){
      if(PROBLEM_TYPE == 1){
        types[i] = GRB_BINARY;
      }
      if(PROBLEM_TYPE == 0){
        types[i] = GRB_CONTINUOUS;
      }

      objective_coeffs[i] = 1.0;
      lower_bounds[i] = 0;
      upper_bounds[i] = 1;
    }

    for(int i=x_start; i<y_start; ++ i){
      if(PROBLEM_TYPE == 1){
        types[i] = GRB_INTEGER;
      }
      if(PROBLEM_TYPE == 0){
        types[i] = GRB_CONTINUOUS;
      }
      objective_coeffs[i] = 0;
      lower_bounds[i] = 0;
      upper_bounds[i] = max_timesteps;
    }

    for(int i=y_start; i<num_variables; ++ i){
      if(PROBLEM_TYPE == 1){
        types[i] = GRB_BINARY;
      }
      if(PROBLEM_TYPE == 0){
        types[i] = GRB_CONTINUOUS;
      }
      objective_coeffs[i] = 0.0;
      lower_bounds[i] = 0;
      upper_bounds[i] = 1.0;
    }

    error = GRBnewmodel(env, &zero_forcing, "forcing", num_variables, objective_coeffs, lower_bounds, upper_bounds, types, NULL);
    if(error){
      cout << "ERROR at: create model" << endl;
      return 1;
    }
    free(objective_coeffs);
    free(lower_bounds);
    free(upper_bounds);
    free(types);

    error = GRBupdatemodel(zero_forcing);
    if(error){
      cout << "ERROR at: first model update" << endl;
      return 1;
    }

    //Now add the constraints

  /*Add a maximal set of disjoint forts*
  {
  bool stopper = false;
  std::set<std::set<int> > forts;
  std::set<int> forbidden_nodes;
  int counter = 0;
  while(!stopper){
    std::vector<double> weights(num_nodes,1.0);
    int size_change = -1*forts.size();
    cout << "forbidden nodes is " << endl;
    for(std::set<int>::iterator jt=forbidden_nodes.begin(); jt!=forbidden_nodes.end(); ++jt){
        cout << *jt << ",";
      }
      cout << endl;
    our_graph.find_forts_LP(env, forts, forbidden_nodes, weights, (double) num_nodes);
    size_change += forts.size();
    if(size_change < 1){
      stopper = true;
    }
  }
    cout << "number of disjoint forts is " << forts.size() << endl;

    //Add fort constraints for start variables
    for(std::set<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
      //mydata.forts.push_back(*it);
      //std::set<int> new_border;
      //bool is_in_border;
      std::vector<int> cind;
      //std::vector<int> cind2;
      std::vector<double> cval;

      for(std::set<int>::iterator jt=it->begin(); jt!=it->end(); ++jt){
        cind.push_back(*jt);
        cval.push_back(-1.0);

      }

      //mydata.fort_borders.push_back(new_border);
      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0, NULL);
      if(error){
        std::cerr << "ERROR: " << error << " disjoint forts constraint" << endl;
        return 1;
      }

    }

  } //*/


    //Constraint 1: each vertex has some method of forcing
    for(int i=0; i<num_nodes; ++i){
      std::vector<int> cind;
      std::vector<double> cval;
      cind.push_back(s_start + i);
      cval.push_back(1.0);
      const int node_degree = our_graph.nodes[i].get_degree();
      for(int j = 0; j < node_degree; ++j){
        cind.push_back(y_start + our_graph.edges.get_edge_number(our_graph.nodes[i].get_adj_list_member(j),i));
        cval.push_back(1.0);
      }


      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_EQUAL, 1.0, NULL);
      if(error){
        cout << "ERROR: " << error << " at constraint 1." << endl;
        return 1;
      }
    }

    //Constraint 2: time constraints for forcing
    for(int i=0; i<num_edges; ++i){
      std::vector<int> cind;
      std::vector<double> cval;
      cind.push_back(y_start + i);
      cval.push_back((double) max_timesteps + 1.0);
      //cind.push_back(y_start +our_graph.edges.get_edge_number(our_graph.edges.get_end2(i), our_graph.edges.get_end1(i)));
      //cval.push_back((double) max_timesteps - 1.0);

      const int edge_end1 = our_graph.edges.get_end1(i);
      const int edge_end2 = our_graph.edges.get_end2(i);

      cind.push_back(x_start + edge_end1);
      cval.push_back(1.0);
      cind.push_back(x_start + edge_end2);
      cval.push_back(-1.0);
      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, (double) max_timesteps, NULL);

      if(error){
        cout << "ERROR at: add loop prevention constraints " << error << endl;
        return 1;
      }
      cind.clear();
      cval.clear();
      cind.push_back(y_start + i);
      cval.push_back((double) max_timesteps + 1.0);
      cind.push_back(x_start + edge_end1);
      cval.push_back(1.0);
      cind.push_back(x_start + edge_end2);
      cval.push_back(-1.0);
      //Now add for 2nd neighbors
      const int neighbor_degree = our_graph.nodes[edge_end1].get_degree();
      for(int j=0; j<neighbor_degree; ++j){
        const int second_neighbor = our_graph.nodes[edge_end1].get_adj_list_member(j);
        if(second_neighbor != edge_end2){
          cind[1] = x_start + second_neighbor;
          error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, (double) max_timesteps, NULL);
          if(error){
            cout << "ERROR at: add loop prevention 2nd neighbor constraints " << error << endl;
            return 1;
          }
        }
      }
    }

  error = GRBupdatemodel(zero_forcing);
    //callback_data* blank_data;
  error = GRBsetcallbackfunc(zero_forcing, ifip::cut_callback, NULL);
  if(error){
	  cout << "ERROR at: GRBsetcallbackfunc" << error << endl;
	  return 1;
	}

  error = GRBsetintparam(GRBgetenv(zero_forcing), "Threads", 1);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "Symmetry", 2);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "MIPFocus", 1);
  error = GRBupdatemodel(zero_forcing);


  error = GRBoptimize(zero_forcing);

  double optimum;
  error = GRBgetdblattr(zero_forcing, GRB_DBL_ATTR_OBJVAL, &optimum);
  // return optimal value
  return nearbyint(optimum);
}