#####################################################
#			Zero Forcing Make File					#
#####################################################
include make.inc

boost_graph:
	g++ -std=gnu++2a cpp/boost_graph.cpp -o boost_graph -I boost_1_75_0
	
zf_test:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o zf_test cpp/zero_forcing.cpp cpp/zf_test.cpp $(CCLNFLAGS)
	
zf_ga:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o zf_ga cpp/zero_forcing.cpp cpp/zf_ga.cpp $(CCLNFLAGS)
	
Source:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o Source cpp/zero_forcing.cpp cpp/Source.cpp $(CCLNFLAGS)
	
wavefront:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o wavefront cpp/zero_forcing.cpp cpp/wavefront.cpp $(CCLNFLAGS)
	
run_wavefront: wavefront
	$(NAUTY)/geng $(SIZE) | ./wavefront
	
run_zf_test: zf_test
	$(NAUTY)/geng $(SIZE) | ./zf_test
	
run_Source: Source
	$(NAUTY)/geng $(SIZE) | ./Source
	
#run_rand_test:
#	$(PYTHON) python/rand_test.py
	
#run_small_test:
#	$(PYTHON) python/small_test.py
	
#zf_test:
#	$(NAUTY)/geng $(SIZE) | $(PYTHON) python/zf_test.py
	
#heuristic:
#	$(NAUTY)/geng $(SIZE) | $(PYTHON) python/heuristic.py
	
py_wavefront:
	$(NAUTY)/geng $(SIZE) | $(PYTHON) python/wavefront.py
	
uninstall:
	@rm -f zf_test
	@rm -f zf_ga
	@rm -f Source
	@rm -f wavefront