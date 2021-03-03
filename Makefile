#####################################################
#			Zero Forcing Make File					#
#####################################################
include make.inc

boost_graph:
	g++ -std=gnu++2a cpp/boost_graph.cpp -o boost_graph -I boost_1_75_0
	
cplex_test:
	g++ -std=gnu++2a cpp/graph6.cpp cpp/cplex_test.cpp -o cplex_test -I $(CPLEX) -I $(CONCERT) -I $(LOCAL)
	
nauty_graph6: cplex_test
	$(NAUTY)/geng $(SIZE) | ./cplex_test
	
run_rand_test:
	$(PYTHON) python/rand_test.py
	
run_small_test:
	$(PYTHON) python/small_test.py
	
uninstall:
	@rm -f boost_graph
	@rm -f graph6
	@rm -f cplex_test