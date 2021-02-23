#####################################################
#			Zero Forcing Make File					#
#####################################################

boost_graph:
	g++ -std=gnu++2a cpp/boost_graph.cpp -o boost_graph -I boost_1_75_0
	
uninstall:
	@rm -f boost_graph