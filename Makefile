#####################################################
#			Zero Forcing Make File					#
#####################################################
include make.inc

boost_graph:
	g++ -std=gnu++2a cpp/boost_graph.cpp -o boost_graph -I boost_1_75_0
	
graph6:
	g++ -std=gnu++2a cpp/graph6.cpp -o graph6
	
nauty_graph6: graph6
	$(NAUTY)/geng $(SIZE) | ./graph6
	
uninstall:
	@rm -f boost_graph
	@rm -f graph6