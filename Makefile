#####################################################
#   Zero-Forcing Make File    						#
#####################################################
include make.inc

install_cubic_test:
	$(CCC) $(CCOPT) $(CCINCFLAGS) $(CCLNFLAGS) -o cubic_test cpp/edg.cpp cpp/wavefront.cpp cpp/genetic_algorithm.cpp cpp/cubic_test.cpp
	
install_nauty_test:
	$(CCC) $(CCOPT) $(CCINCFLAGS) $(CCLNFLAGS) -o nauty_test cpp/graph6.cpp cpp/wavefront.cpp cpp/genetic_algorithm.cpp cpp/nauty_test.cpp
	
install_watts_strogatz_test:
	$(CCC) $(CCOPT) $(CCINCFLAGS) $(CCLNFLAGS) -o watts_strogatz_test cpp/edg.cpp cpp/wavefront.cpp cpp/genetic_algorithm.cpp cpp/watts_strogatz_test.cpp
	
install_star_test:
	$(CCC) $(CCOPT) $(CCINCFLAGS) $(CCLNFLAGS) -o star_test cpp/edg.cpp cpp/infection_ip.cpp cpp/genetic_algorithm.cpp cpp/star_test.cpp
	
run_all: install_cubic_test install_watts_strogatz_test install_star_test
	./cubic_test
	./watts_strogatz_test
	./star_test
	
run_nauty_all: install_nauty_test
	$(NAUTY)/geng $(ORDER) | ./nauty_test
	
run_nauty_rand: install_nauty_test
	$(NAUTY)/genrang $(ORDER) 100 -g -P1/2 | ./nauty_test

uninstall:
	@rm -f cubic_test
	@rm -f watts_strogatz_test
	@rm -f star_test