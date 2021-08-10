#include <iostream>
#include <fstream>
#include <ip.hpp>
#include <ext_cover_ip.hpp>
#define STOPPING_TIME 7200
#define PROBLEM_TYPE GRB_BINARY
using namespace std;
using namespace ip;
using namespace ecip;
/* cut_callback */
int __stdcall ecip::cut_callback(GRBmodel *model, void *cbdata, int where, void *usrdata){

  if(where != GRB_CB_POLLING){
    struct callback_data* mydata = (struct callback_data *) usrdata;
    	auto current = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> running_time = current - mydata->chrono_start_time;

    /*clock_t current_time = clock();
    double running_time = (current_time -mydata->start_time)/static_cast<double>(CLOCKS_PER_SEC);

      if(running_time > STOPPING_TIME){
	      GRBterminate(model);
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
    /*for(int i=0; i<mydata->num_variables; ++i){
      weights[i] = incumbent_solution[i];
      }*/

    //Find the end nodes of the current incumbent
    ip::node_subset new_subset;
    for(int i=0; i<num_nodes; ++i){
      if(incumbent_solution[i] > 1 - FLOAT_TOL){
        new_subset.add_to_subset(i);
      }
    }
    for(int j=0; j<mydata->node_subsets.size(); ++j){
      if(incumbent_solution[mydata->a_start+j] > 1 - FLOAT_TOL){
        for(std::set<int>::iterator it = mydata->node_subsets[j].get_subset_begin(); it != mydata->node_subsets[j].get_subset_end(); ++it){
          new_subset.add_to_subset(*it);
        }
        /*new_subset.add_to_subset(j);
        const int num_neighbors = mydata->our_graph.nodes[j].get_degree();
        for(int i=0; i <num_neighbors-1; ++i){
          new_subset.add_to_subset(mydata->our_graph.nodes[j].get_adj_list_member(i));
        }*/
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
    auto time_start = std::chrono::high_resolution_clock::now();

    std::set<std::set<int> > forts;
    // mydata->our_graph.find_fort_closure(forts, weights);
    //mydata->our_graph.find_fort_closure(forts, end_weights);
    //mydata->our_graph.find_fort_greedy_style_2(forts, weights);
    //mydata->our_graph.find_fort_greedy_style_2(forts, end_weights);

    //mydata->our_graph.find_multiple_minimal_border_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);
    std::set<int> forbidden_nodes;

    for(std::set<int>::iterator it = new_subset.get_closure_begin(); it != new_subset.get_closure_end(); ++it){
      forbidden_nodes.insert(*it);
    }
    /*//std::set<int> new_fort;
    //for(int i=0; i<mydata->num_nodes; ++i){
    //if(!new_subset.is_in_closure(i)){
	//cout << "inserting " << i << endl;
    //	new_fort.insert(i);
    //}
    //}
    //cout << "inserting fort" << endl;
    //if(new_fort.size() > 0){
    //forts.insert(new_fort);
    //}
    //mydata->our_graph.find_minimal_fort_LP(mydata->env, forts, forbidden_nodes, weights, 1.0);
    //Get forbidden nodes as intersection of previous forts */

    //std::vector<double> one_weights(num_nodes, 1.0);
    mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, forbidden_nodes, weights, (double) num_nodes, mydata->num_forts_variable_in); //Current best
    if(forts.size() > 0){
      //cout << "FOUND A NEW FORT" << endl;
      //mydata->num_forts++;
    }
    //mydata->our_graph.find_minimal_border_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0); //Current best

    //mydata->our_graph.find_max_intersection_fort_LP(mydata->env, forts, mydata->num_forts_variable_in, fort_size, weights, 1.0);
    //mydata->our_graph.find_maximal_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);
    //mydata->our_graph.find_minimal_fort_LP(mydata->env, forts, mydata->nodes_in_cover, weights, 1.0);
    //mydata->our_graph.find_forts_LP(mydata->env, forts, mydata->nodes_in_cover, fort_number_weights, 1.0);

    //Add the fort constraints for start vertices
    for(std::set<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
    std::set<int> A;
    std::set<int> B;
    int v = -1;
    //int v = mydata->our_graph.determine_facet_LP(mydata->env, *it, A, B);
      if(v >= 0){
        std::vector<int> cind;
        std::vector<double> cval;
        cind.push_back(v);
        cval.push_back(1.0);
        std::set<std::pair<int,double> > a_variables;

        for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
          cind.push_back(*jt);
          cval.push_back(1.0);

          for(int j=0; j<mydata->node_subsets.size(); ++j){
          if(mydata->node_subsets[j].is_in_closure(*jt)){
            std::pair<int, double> new_pair;
            new_pair.first = j;
            if(A.count(*jt) == 1 && B.count(*jt) == 1){
              new_pair.second = 2.0;
            }
            else{
              new_pair.second = 1.0;
            }
            mydata->num_forts_variable_in[num_nodes + j]++;
          }
          }


        }
        for(std::set<std::pair<int,double> >::iterator kt = a_variables.begin(); kt != a_variables.end(); ++kt){
          cind.push_back(mydata->a_start + kt->first);
          cval.push_back(kt->second);
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
      std::vector<int> cind;
      //std::vector<int> cind2;
      std::vector<double> cval;
      std::set<int> a_variables;
      //cout << "adding constraint with: ";
      for(std::set<int>::iterator jt = it->begin(); jt != it->end(); ++jt){
        //mydata->forts_node_in[*jt].push_back(mydata->forts.size()-1);
        for(int i=0; i<mydata->our_graph.nodes[*jt].get_degree(); ++i){
          if(it->count(mydata->our_graph.nodes[*jt].get_adj_list_member(i)) != 1){
            //cout << "node is on border " << endl;
            //new_border.insert(*jt);
            mydata->num_forts_variable_in[*jt]++;
            cind.push_back(*jt);
            //cind2.push_back(mydata->end_variable_start + *jt);
            cval.push_back(-1.0);
            break;
          }
        }
        for(int j=0; j<mydata->node_subsets.size(); ++j){
          if(mydata->node_subsets[j].is_in_closure(*jt)){
            a_variables.insert(j);
            mydata->num_forts_variable_in[num_nodes + j]++;
          }
        }

      }
      for(std::set<int>::iterator kt = a_variables.begin(); kt != a_variables.end(); ++kt){
        cind.push_back(mydata->a_start + *kt);
        cval.push_back(-1.0);
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
      error = GRBcblazy(cbdata, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0);
      if(error){
	cout << "ERROR at: add callback fort constraint " << error << endl;
	return 1;
      }
      mydata->num_forts++;

    }
    }

	auto time_end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> running_time = time_end - time_start;
	mydata->time_in_LP += running_time.count();

    /*for(int i=0;i<2*num_nodes; ++i){
      error = GRBsetintattrelement(model, GRB_INT_ATTR_BRANCHPRIORITY, i, -1*mydata->num_forts_variable_in[i]);
      if(error){
      	cout << "ERROR at: changing branchpriority " << error << endl;
        return 1;
      }
    }*/
  }


  return 0;
};
/* zero_forcing */
int ecip::zero_forcing(ip::graph &our_graph)
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

    std::vector<node_subset> node_subsets;
    std::vector<int> subset_hashes;
    for(int i=0; i<num_nodes; ++i){
      node_subset new_subset;
      new_subset.add_to_subset(i);
      const int num_neighbors = our_graph.nodes[i].get_degree();
      for(int j=0; j<num_neighbors-1; ++j){
        new_subset.add_to_subset(our_graph.nodes[i].get_adj_list_member(j));
      }
      new_subset.set_weight();
      new_subset.update_closure(our_graph.nodes);
      const int new_hash = new_subset.create_hash();
      bool same_closure = false;
      for(int j=0; j < node_subsets.size(); ++j){
        //cout << "Hash for " << j << " is " << node_subsets[j].create_hash() << endl;
        if(new_hash == subset_hashes[j]){
          //cout << "CLOSURES ARE THE SAME" << endl;
          same_closure = true;
          break;
        }
      }
      if(!same_closure){
        node_subsets.push_back(new_subset);
        subset_hashes.push_back(new_hash);
      }
    }
    /*Create the model*******************************************************/
    //Start with number of forcing chains = to heuristic branch_width value
    //Create the model with 3N*branch_width binary variables
    const int num_variables = num_nodes + node_subsets.size();
    const int a_start = num_nodes;

  //Set Callback function for adding cuts
    callback_data mydata;
    mydata.previous_size = 0;
  mydata.counter = 0;
    mydata.a_start = a_start;
  mydata.min_degree = min_degree;
    mydata.node_subsets = node_subsets;
  //mydata.edges = edges;
  //mydata.nodes = nodes;
    mydata.num_forts = 0;
  mydata.num_forts_variable_in = std::vector<int>(num_variables, 0);
  mydata.num_nodes = num_nodes;
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
    for(int i=0; i<a_start; ++ i){
      types[i] = PROBLEM_TYPE;
      //types[i] = GRB_CONTINUOUS;
      //Objective sums all nodes at first time stage
      objective_coeffs[i] = 1;
      lower_bounds[i] = 0;
      upper_bounds[i] = 1;
    }
    for(int i=a_start; i<num_variables; ++ i){
      types[i] = PROBLEM_TYPE;
      //types[i] = GRB_CONTINUOUS;
      //Objective sums all nodes at first time stage
      objective_coeffs[i] = (double) node_subsets[i-a_start].get_weight();
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

    /*Constraint 1: Add fort constraints that do not contain subset variables*/
  {
  std::set<std::set<int> > forts;
  for(int i=0; i<node_subsets.size(); ++i){
      std::set<int> new_fort;
      //std::set<int> forbidden_nodes;
      //std::vector<double> weights(num_nodes, 0);
      //std::vector<int> blank;
      for(int j=0; j<num_nodes; ++j){
        if(!node_subsets[i].is_in_closure(j)){
          //forbidden_nodes.insert(i);
          new_fort.insert(j);
        }
      }
      //our_graph.find_minimal_border_fort_LP(env, forts, forbidden_nodes, weights, (double) num_nodes,  blank);
      if(new_fort.size() >=1 ){
	forts.insert(new_fort);
      }
      }


    for(std::set<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
      mydata.num_forts++;
      //mydata.forts.push_back(*it);
      //std::set<int> new_border;
      std::vector<int> cind;
      //std::vector<int> cind2;
      std::vector<double> cval;
      std::set<int> a_variables;
      for(std::set<int>::iterator jt=it->begin(); jt!=it->end(); ++jt){
        for(int i=0; i<our_graph.nodes[*jt].get_degree(); ++i){
          if(it->count(our_graph.nodes[*jt].get_adj_list_member(i)) != 1){
            //cout << "node is on border " << endl;
            //new_border.insert(*jt);
            //mydata.num_forts_variable_in[*jt]++;
            cind.push_back(*jt);
            //cind2.push_back(end_variable_start+*jt);
            cval.push_back(-1.0);
            break;
          }
        }

        for(int j=0; j<mydata.node_subsets.size(); ++j){
          if(mydata.node_subsets[j].is_in_closure(*jt)){
            //cout << "node is in closure " << endl;
            a_variables.insert(j);
          }
        }
      }

      for(std::set<int>::iterator kt = a_variables.begin(); kt != a_variables.end(); ++kt){
        cind.push_back(mydata.a_start + *kt);
        //cind2.push_back(ae_start + *kt);
        //mydata.num_forts_variable_in[mydata.a_start + *kt]++;
        //cind2.push_back(ae_start + *kt);
        cval.push_back(-1.0);
      }

      //mydata.fort_borders.push_back(new_border);
      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0, NULL);
      if(error){
        std::cerr << "ERROR: " << error << " disjoint forts constraint" << endl;
        return 1;
      }

    }
  }
  //*/

   /*Constraint 3: a variables are 0 if they are part of another positive a variable *
  for(int i=0; i<num_nodes; ++i){
    std::vector<int> cind;
    std::vector<double> cval;
    for(int j=0; j<num_nodes; ++j){
      for(std::set<int>::iterator jt = node_subsets[j].get_subset_begin(); jt != node_subsets[j].get_subset_end(); ++jt){
      if(node_subsets[i].is_in_closure(*jt)){
        cind.push_back(a_start+j);
        cval.push_back(1.0);
        break;
      }
    }
    error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 1, NULL);
      if(error){
        std::cerr << "ERROR: " << error << "at a variable disjointness constraint" << endl;
        return 1;
      }
  }//*/



      /* start variables are only nonzero if they are outside of a variable closures*
    for(int i=0; i<num_nodes; ++i){
      for(int j=0; j<num_nodes; ++j){
          if(node_subsets[j].is_in_closure(i)){
                          std::vector<int> cind;
          std::vector<double> cval;
          cind.push_back(i);
          cval.push_back(1.0);
            cind.push_back(a_start+j);
            cval.push_back(1.0);
                            error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, 1, NULL);
          if(error){
            std::cerr << "ERROR: " << error << "no start in closures" << endl;
            return 1;
          }
          }
        }
      }

    }


         /*Constraint 2: There must be one positive a variable*****************
    {
      std::vector<int> cind;
      std::vector<double> cval;
      for(int i=0; i<num_nodes; ++i){
        cind.push_back(a_start+i);
        cval.push_back(1.0);
      }
      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_GREATER_EQUAL, 1, NULL);
      if(error){
	std::cerr << "ERROR: " << error << "at must be one positive a variable" << endl;
	return 1;
      }
    }

    /*Add a maximal set of disjoint forts*/
  {
  bool stopper = false;
  std::set<std::set<int> > forts;
  std::set<int> forbidden_nodes;
  int counter = 0;
  while(!stopper){
    std::vector<double> weights(num_nodes,1.0);
    int size_change = -1*forts.size();
    our_graph.find_forts_LP(env, forts, forbidden_nodes, weights, (double) num_nodes);
    size_change += forts.size();
    if(size_change < 1){
      stopper = true;
    }
  }

    //Add fort constraints for start variables
    for(std::set<std::set<int> >::iterator it = forts.begin(); it != forts.end(); ++it){
      //mydata.forts.push_back(*it);
      //std::set<int> new_border;
      std::vector<int> cind;
      //std::vector<int> cind2;
      std::vector<double> cval;
      std::set<int> a_variables;
      for(std::set<int>::iterator jt=it->begin(); jt!=it->end(); ++jt){
        for(int i=0; i<our_graph.nodes[*jt].get_degree(); ++i){
          if(it->count(our_graph.nodes[*jt].get_adj_list_member(i)) != 1){
            //cout << "node is on border " << endl;
            //new_border.insert(*jt);
            //mydata.num_forts_variable_in[*jt]++;
            cind.push_back(*jt);
            //cind2.push_back(end_variable_start+*jt);
            cval.push_back(-1.0);
            break;
          }
        }
        for(int j=0; j<mydata.node_subsets.size(); ++j){
          if(mydata.node_subsets[j].is_in_closure(*jt)){
            a_variables.insert(j);
          }
        }

      }
      for(std::set<int>::iterator kt = a_variables.begin(); kt != a_variables.end(); ++kt){
        cind.push_back(mydata.a_start + *kt);
        cval.push_back(-1.0);
      }
      //mydata.fort_borders.push_back(new_border);
      error = GRBaddconstr(zero_forcing, cind.size(), &cind[0], &cval[0], GRB_LESS_EQUAL, -1.0, NULL);
      if(error){
        std::cerr << "ERROR: " << error << " disjoint forts constraint" << endl;
        return 1;
      }
      mydata.num_forts++;

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

  error = GRBsetcallbackfunc(zero_forcing, ecip::cut_callback, (void *) &mydata);
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