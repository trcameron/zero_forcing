#####################################################
#			Zero Forcing Make File					#
#####################################################
include make.inc

boost_graph:
	g++ -std=gnu++2a cpp/boost_graph.cpp -o boost_graph -I boost_1_75_0
	
heuristic:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o heuristic cpp/zero_forcing.cpp cpp/heuristic.cpp $(CCLNFLAGS)
	
run_heuristic: heuristic
	@echo "Min Ratio, Avg Ratio, Max Ratio"
	@$(NAUTY)/geng $(N) -q | ./heuristic

zf_test_corr:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o zf_test_corr cpp/zero_forcing.cpp cpp/zf_test_corr.cpp $(CCLNFLAGS)
	
zf_test_time:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o zf_test_time cpp/zero_forcing.cpp cpp/zf_test_time.cpp $(CCLNFLAGS)
	
all_zf_test_corr: zf_test_corr
	@echo "Wave % Correct, GA % Correct"
	@$(NAUTY)/geng $(N) -q | ./zf_test_corr
	
rand_zf_test_corr: zf_test_corr
	@echo "Wave % Correct, GA % Correct"
	@$(NAUTY)/genrang $(N) $(NUM) -g -P20/100 | ./zf_test_corr
	@$(NAUTY)/genrang $(N) $(NUM) -g -P30/100 | ./zf_test_corr
	@$(NAUTY)/genrang $(N) $(NUM) -g -P40/100 | ./zf_test_corr
	@$(NAUTY)/genrang $(N) $(NUM) -g -P50/100 | ./zf_test_corr
	@$(NAUTY)/genrang $(N) $(NUM) -g -P60/100 | ./zf_test_corr
	@$(NAUTY)/genrang $(N) $(NUM) -g -P70/100 | ./zf_test_corr
	@$(NAUTY)/genrang $(N) $(NUM) -g -P80/100 | ./zf_test_corr
	
all_zf_test_time: zf_test_time
	@echo "Avg IP Time, Avg Wave Time, Avg GA Time"
	@$(NAUTY)/geng $(N) -q | ./zf_test_time
	
rand_zf_test_time: zf_test_time
	@echo "Avg IP Time, Avg Wave Time, Avg GA Time"
	@$(NAUTY)/genrang $(N) $(NUM) -g -P20/100 | ./zf_test_time
	@$(NAUTY)/genrang $(N) $(NUM) -g -P30/100 | ./zf_test_time
	@$(NAUTY)/genrang $(N) $(NUM) -g -P40/100 | ./zf_test_time
	@$(NAUTY)/genrang $(N) $(NUM) -g -P50/100 | ./zf_test_time
	@$(NAUTY)/genrang $(N) $(NUM) -g -P60/100 | ./zf_test_time
	@$(NAUTY)/genrang $(N) $(NUM) -g -P70/100 | ./zf_test_time
	@$(NAUTY)/genrang $(N) $(NUM) -g -P80/100 | ./zf_test_time
	
#run_rand_test:
#	$(PYTHON) python/rand_test.py
	
#run_small_test:
#	$(PYTHON) python/small_test.py
	
#zf_test:
#	$(NAUTY)/geng $(SIZE) | $(PYTHON) python/zf_test.py
	
#heuristic:
#	$(NAUTY)/geng $(SIZE) | $(PYTHON) python/heuristic.py
	
py_wavefront:
	$(NAUTY)/geng $(N) | $(PYTHON) python/wavefront.py
	
uninstall:
	@rm -f zf_test_corr
	@rm -f zf_test_time
	@rm -f heuristic