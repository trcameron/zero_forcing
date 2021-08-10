#include <iostream>
#include <fstream>
#include <ip.hpp>
#include <max_closure_ip.hpp>
#define STOPPING_TIME 7200
using namespace std;
using namespace ip;
using namespace mcip;
/* cut_callback */
int __stdcall mcip::cut_callback(GRBmodel *model, void *cbdata, int where, void *usrdata){

  if(where != GRB_CB_POLLING){
        struct callback_data* mydata = (struct callback_data *) usrdata;
	//clock_t current_time = clock();
	auto current = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> running_time = current - mydata->chrono_start_time;
	//double running_time = (current_time - mydata->start_time)/static_cast<double>(CLOCKS_PER_SEC) ;
    /*int error;
      error = GRBcbget(cbdata, where, GRB_CB_RUNTIME, &running_time);

      if(error){
        cerr << "ERROR: " << error << " getting runtime." << endl;
	}*/
	if(running_time > std::chrono::seconds(STOPPING_TIME)){
	      GRBterminate(model);
      }
  }

  if(where == GRB_CB_MIPSOL){
    struct callback_data* mydata = (struct callback_data *) usrdata;
    int error=0;

    //cout << "found new incumbent" << endl;
    const int num_nodes = mydata->num_nodes;
    int fort_size = mydata->forts.size();
    //cout << "forts.size() is " << fort_size << endl;

    double* incumbent_solution = (double*) calloc(mydata->num_variables, sizeof(double));
    error = GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, incumbent_solution);
    std::vector<double> weights(mydata->num_variables, 0);
    for(int i=0; i<mydata->num_variables; ++i){
      weights[i] = incumbent_solution[i];
    }

    //Find the end nodes of the current incumbent
    ip::node_subset new_subset;
    for(int i=0; i<num_nodes; ++i){
      if(incumbent_solution[i] > 1 - FLOAT_TOL){
        new_subset.add_to_subset(i);
      }
    }
    for(int j=0; j<mydata->node_subsets.size(); ++j){
      if(incumbent_solution[mydata->a_start+j] > 1 - FLOAT_TOL){
        new_subset.add_to_subset(j);
        const int num_neighbors = mydata->our_graph.nodes[j].get_degree();
        for(int i=0; i <num_neighbors-1; ++i){
          new_subset.add_to_subset(mydata->our_graph.nodes[j].get_adj_list_member(i));
        }
      }
    }


    std::set<int> end_nodes;
    new_subset.set_weight();
    new_subset.update_closure_get_ends(mydata->our_graph.nodes, end_nodes);
    /*std::vector<double> end_weights(num_nodes, 0);
    for(int i=0; i<num_nodes; ++i){
      if(end_nodes.count(i) == 1){
        end_weights[i] = 1.0;
      }
    }*/


    //Find the violated fort constraints


    /*std::vector<double> fort_number_weights(mydata->num_variables, 0);
    for(int i=0; i<mydata->num_variables; ++i){
      fort_number_weights[i]  = weights[i] + static_cast<double>(mydata->num_forts_variable_in[i])/(100.0*fort_size);
    }*/
    //clock_t time_start = clock();
    auto time_start = std::chrono::high_resolution_clock::now();
    std::set<std::set<int> > forts;
    //mydata->our_graph.find_fort_closure(forts, weights);
    //mydata->our_graph.find_fort_closure(forts, end_weights);
    //mydata->our_graph.find_fort_greedy_style_2(forts, weights);
    //mydata->our_graph.find_fort_greedy_style_2(forts, end_weights);

    //mydata->our_graph.find_multiple_minimal_border_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);
    std::set<int> forbidden_nodes;

    for(std::set<int>::iterator it = new_subset.get_closure_begin(); it != new_subset.get_closure_end(); ++it){
      forbidden_nodes.insert(*it);
    }
    mydata->our_graph.find_fort_max_complement(forts, forbidden_nodes);
       //Get forbidden nodes as intersection of previous forts

    //std::vector<double> one_weights(num_nodes, 1.0);
    //mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, forbidden_nodes, weights, 1.0, mydata->num_forts_variable_in); //Current best

    //mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0); //Current best

    //mydata->our_graph.find_max_intersection_fort_LP(mydata->env, forts, mydata->num_forts_variable_in, fort_size, weights, 1.0);
    //mydata->our_graph.find_maximal_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);
    //mydata->our_graph.find_minimal_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);
    //mydata->our_graph.find_forts_LP(mydata->env, forts, mydata->nodes_in_cover, fort_number_weights, 1.0);

    //Add the fort constraints for start vertices
    for(std::set<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
		int v = mydata->our_graph.determine_facet_AB_LP(mydata->env, *it);

      //int v = -1;
      if(v >= 0){
        std::vector<int> cind;
        std::vector<double> cval;
        cind.push_back(v);
        cval.push_back(1.0);
        for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
          cind.push_back(*jt);
          cval.push_back(1.0);
        }
        error = GRBcblazy(cbdata, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, 2);
        if(error){
          cout << "ERROR at: add callback gomory constraint " << error << endl;
          return 1;
        }
      }
      else{
      //mydata->forts.push_back(*it);
      //std::set<int> new_border;
      bool is_in_border;
      std::vector<int> cind;
      //std::vector<int> cind2;
      std::vector<double> cval;
      std::set<int> a_variables;
      //cout << "adding constraint with: ";
      int min_degree = INT_MAX;
      for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
        const int degree = mydata->our_graph.nodes[*jt].get_degree();
        if(degree < min_degree){
          min_degree = degree;
        }
      }


      for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
        //mydata->forts_node_in[*jt].push_back(mydata->forts.size()-1);
        is_in_border = false;
        for(int i=0; i<mydata->our_graph.nodes[*jt].get_degree(); ++i){
          if(it->count(mydata->our_graph.nodes[*jt].get_adj_list_member(i)) != 1){
            //cout << "node is on border " << endl;
            //new_border.insert(*jt);
            is_in_border = true;
            //mydata->num_forts_variable_in[*jt]++;
            //cind.push_back(get_x_v_S_number(num_nodes,*jt));
            //cind2.push_back(mydata->end_variable_start + *jt);
            //cval.push_back(-1.0);
            break;
          }
        }
        if(is_in_border){
          cind.push_back(*jt);
          cval.push_back(-1.0);
        }
        else{
          cind.push_back(*jt);
          cval.push_back(-1.0/min_degree);
        }
      }


        //cind2.push_back(mydata->end_variable_start + *jt);


	//Find out the minimum "degree" of the nodes adjacent to each node in the fort
	//fort_degree is the number of neighbors inside the fort
	//min_out_degree is the number of neighbors outside the fort
	/*
	const int node_degree = mydata->our_graph.nodes[*jt].get_degree();
	int fort_degree = 0;
	int min_out_degree = num_nodes;
	for(int i=0; i<node_degree; ++i){
	  const int neighbor = mydata->our_graph.nodes[*jt].get_adj_list_member(i);
	  //If not in fort, check neighbors degree of nodes in fort
	  if(it->count(neighbor) == 0){
	    const int neighbor_degree = mydata->our_graph.nodes[neighbor].get_degree();
	    int degree = 0;
	    for(int j=0; j<neighbor_degree; ++j){
	      if(it->count(mydata->our_graph.nodes[neighbor].get_adj_list_member(j))==1){
          degree++;
	      }
	    }
	    if(degree < min_out_degree){
	      min_out_degree = degree;
	    }
	  }
	  else{
	    fort_degree++;
	  }
	}
	if(fort_degree == 0){
	  fort_degree++;
	}
	if(min_out_degree < 2){
	  std::cerr << "ERROR: Fort is not actually a fort!" << endl;
	}
        for(int i=0; i<it->size(); ++i){
         // cval.push_back((fort_degree > (min_out_degree-1)) ? -1.0/((double)min_out_degree-1) : -1.0/(double)fort_degree);
        }
  */
	//cout << endl;
	//cout << "cind.size() is " << cind.size() << endl;

      //mydata->fort_borders.push_back(new_border);
      mydata->num_forts++;
      error = GRBcblazy(cbdata, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0);
      if(error){
	cout << "ERROR at: add callback fort constraint " << error << endl;
	return 1;
      }
      }
    }

        //clock_t time_end = clock();
        auto time_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> running_time = time_end - time_start;
    mydata->time_in_LP += running_time.count();
  }

  return 0;
};
/* zero_forcing */
int mcip::zero_forcing(ip::graph &our_graph)
{
    GRBenv *env = NULL;
    GRBmodel *zero_forcing = NULL;
    int error = 0;
    int num_nodes = our_graph.nodes.size();
    int num_edges = our_graph.edges.get_size();

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
    //Start with number of forcing chains = to heuristic branch_width value
    //Create the model with 3N*branch_width binary variables
    const int num_variables = num_nodes;


  //Set Callback function for adding cuts
    callback_data mydata;
    mydata.previous_size = 0;
  mydata.counter = 0;
  mydata.min_degree = min_degree;
  //mydata.edges = edges;
  //mydata.nodes = nodes;
  mydata.num_forts_variable_in = std::vector<int>(num_variables, 0);
  mydata.num_nodes = num_nodes;
  mydata.num_forts = 0;
  mydata.num_variables = num_variables;
  //mydata.our_graph = our_graph;
  mydata.env = env;
    mydata.our_graph.nodes = our_graph.nodes;
    const int edge_list_size = our_graph.edges.get_size();
    for(int i=0; i<edge_list_size; ++i){
      mydata.our_graph.edges.add_edge(our_graph.edges.get_end1(i), our_graph.edges.get_end2(i));
    }
  mydata.forts_node_in.resize(num_nodes);
    char* types = (char*) calloc(num_variables, sizeof(char));
    double* objective_coeffs = (double*) calloc(num_variables, sizeof(double));
    double* upper_bounds = (double*) calloc(num_variables, sizeof(double));
    double* lower_bounds = (double*) calloc(num_variables, sizeof(double));

    //Variables for first time stage (x_v,0 variables)
    for(int i=0; i<num_variables; ++ i){
      types[i] = GRB_BINARY;
      //types[i] = GRB_CONTINUOUS;
      //Objective sums all nodes at first time stage
      objective_coeffs[i] = 1;
      lower_bounds[i] = 0;
      upper_bounds[i] = 1;

    }


    /***********************/
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


  /*add constraints for unnessential pairs*
  std::vector<std::pair<int,int> > unnessential_pairs;
  unnessential_pairs.push_back(std::make_pair(40,12));
  unnessential_pairs.push_back(std::make_pair(20,18));
  unnessential_pairs.push_back(std::make_pair(6,30));
  unnessential_pairs.push_back(std::make_pair(29,31));
  unnessential_pairs.push_back(std::make_pair(21,44));
  unnessential_pairs.push_back(std::make_pair(45,46));
  unnessential_pairs.push_back(std::make_pair(10,23));
  unnessential_pairs.push_back(std::make_pair(43,8));
  unnessential_pairs.push_back(std::make_pair(17,39));
  unnessential_pairs.push_back(std::make_pair(4,36));
  unnessential_pairs.push_back(std::make_pair(7,48));
  unnessential_pairs.push_back(std::make_pair(56,55));
  unnessential_pairs.push_back(std::make_pair(53,54));
  unnessential_pairs.push_back(std::make_pair(27,28));
  unnessential_pairs.push_back(std::make_pair(42,0));
  unnessential_pairs.push_back(std::make_pair(57,58));
  unnessential_pairs.push_back(std::make_pair(59,14));
  unnessential_pairs.push_back(std::make_pair(15,38));
  unnessential_pairs.push_back(std::make_pair(16,3));
  unnessential_pairs.push_back(std::make_pair(34,37));
  unnessential_pairs.push_back(std::make_pair(24,35));

  for(std::vector<std::pair<int,int> >::iterator it = unnessential_pairs.begin(); it != unnessential_pairs.end(); ++it){
    std::vector<int> cind;
    std::vector<double> cval;
    cind.push_back(it->first);
    cval.push_back(1.0);
    cind.push_back(it->second);
    cval.push_back(1.0);

    error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 1.0, NULL);
    if(error){
      std::cerr << "ERROR: " << error << " at unnessential pairs constraint" << endl;
      return 1;
    }

  }//*/

      /*Add a maximal set of disjoint forts*/
  {
  bool stopper = false;
  std::set<std::set<int> > forts;
  std::set<int> forbidden_nodes;
  int counter = 0;
  while(!stopper){
    std::vector<double> weights(num_nodes,1.0);
    int size_change = -1*forts.size();
      //our_graph.find_forts_LP(env, forts, forbidden_nodes, weights, (double) num_nodes);
      our_graph.find_fort_max_complement(forts, forbidden_nodes);
      for(std::set<std::set<int> >::iterator st = forts.begin(); st != forts.end(); ++st){
        for(std::set<int>::iterator it = st->begin(); it != st->end(); ++it){
	  forbidden_nodes.insert(*it);
        }
      }
    size_change += forts.size();
    if(size_change < 1){
      stopper = true;
    }
  }

    //Add fort constraints for start variables
    for(std::set<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
      //mydata.forts.push_back(*it);
      //std::set<int> new_border;
      bool is_in_border;
      std::vector<int> cind;
      //std::vector<int> cind2;
      std::vector<double> cval;
      mydata.num_forts++;
      /*int min_degree = INT_MAX;
      for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
        const int degree = our_graph.nodes[*jt].get_degree();
        cout << "degree is " << degree << endl;
        if(degree < min_degree){
          min_degree = degree;
          cout << "min degree is " << min_degree << endl;
        }
      }*/

      for(std::set<int>::iterator jt=it->begin(); jt!=it->end(); ++jt){
        cind.push_back(*jt);
        cval.push_back(-1.0);
        /*is_in_border = false;
        for(int i=0; i<our_graph.nodes[*jt].get_degree(); ++i){
          if(it->count(our_graph.nodes[*jt].get_adj_list_member(i)) != 1){
            //cout << "node is on border " << endl;
            //new_border.insert(*jt);
            is_in_border = true;
            //mydata.num_forts_variable_in[*jt]++;
            //cind.push_back(get_x_v_S_number(num_nodes,*jt));
            //cind2.push_back(end_variable_start+*jt);
            //cval.push_back(-1.0);
            break;
          }
        }
        if(is_in_border){
          cind.push_back(*jt);
          cval.push_back(-1.0);
        }
        else{
          cind.push_back(*jt);
          cval.push_back(-1.0/min_degree);
        }*/
      }

      //mydata.fort_borders.push_back(new_border);
      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0, NULL);
      if(error){
        std::cerr << "ERROR: " << error << " disjoint forts constraint" << endl;
        return 1;
      }

    }

  } //*/
//*/

 //*/
	//End Gurobi Model***************************************************************************
  //error = GRBwrite(zero_forcing, "zero_forcing.lp");
	//Solve the model
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "PreSolve", 0);
  error = GRBsetintparam(GRBgetenv(zero_forcing), "PreCrush", 1);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "Method", 1);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "MIPFocus", 1);
 // error = GRBsetintparam(GRBgetenv(zero_forcing), "RINS", 5);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "BranchDir", 1);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "Cuts", 0);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "GomoryPasses", 5);
  error = GRBsetintparam(GRBgetenv(zero_forcing), "LazyConstraints", 1);
  error = GRBsetintparam(GRBgetenv(zero_forcing), "Threads", 1);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "Symmetry", 2);
  //error = GRBsetintparam(GRBgetenv(zero_forcing), "OutputFlag", 0);


  //mydata.orbits = orbits;
  error = GRBsetcallbackfunc(zero_forcing, cut_callback, (void *) &mydata);
  if(error){
	  cout << "ERROR at: GRBsetcallbackfunc" << error << endl;
	  return 1;
	}

  error = GRBupdatemodel(zero_forcing);
  mydata.chrono_start_time = std::chrono::high_resolution_clock::now();
  error = GRBoptimize(zero_forcing);
  double optimum;
  error = GRBgetdblattr(zero_forcing, GRB_DBL_ATTR_OBJVAL, &optimum);

  error = GRBfreemodel(zero_forcing);
    if(error){
    cout << "ERROR at: free model" << endl;
    return 1;
  }
  GRBfreeenv(env);
  // return optimal value
  return nearbyint(optimum);
}