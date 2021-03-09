#####################################################
#			Zero Forcing Make File					#
#####################################################
include make.inc

boost_graph:
	g++ -std=gnu++2a cpp/boost_graph.cpp -o boost_graph -I boost_1_75_0
	
cplex_test:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o cplex_test cpp/graph6.cpp cpp/cplex_test.cpp $(CCLNFLAGS)
	#g++ -std=gnu++2a -m64 -fPIC cpp/graph6.cpp cpp/cplex_test.cpp -o cplex_test -I $(CPLEXI) -I $(CONCERTI) -I $(LOCAL) -L $(CPLEXL) -L $(CONCERTL)
	
nauty_graph6: cplex_test
	$(NAUTY)/geng $(SIZE) | ./cplex_test
	
run_rand_test:
	$(PYTHON) python/rand_test.py
	
run_small_test:
	$(PYTHON) python/small_test.py
	
uninstall:
	@rm -f cplex_test