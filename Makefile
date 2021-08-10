#####################################################
#   Zero-Forcing Make File    						#
#####################################################
include make.inc

install_cubic_test:
	$(CCC) $(CCOPT) $(CCINCFLAGS) $(CCLNFLAGS) -o cubic_test cpp/edg.cpp cpp/wavefront.cpp cpp/genetic_algorithm.cpp cpp/cubic_test.cpp
	
install_watts_strogatz_test:
	$(CCC) $(CCOPT) $(CCINCFLAGS) $(CCLNFLAGS) -o watts_strogatz_test cpp/edg.cpp cpp/wavefront.cpp cpp/genetic_algorithm.cpp cpp/watts_strogatz_test.cpp
	
install_star_test:
	$(CCC) $(CCOPT) $(CCINCFLAGS) $(CCLNFLAGS) -o star_test cpp/edg.cpp cpp/infection_ip.cpp cpp/genetic_algorithm.cpp cpp/star_test.cpp
	
run_all: install_cubic_test install_watts_strogatz_test install_star_test
	./cubic_test
	./watts_strogatz_test
	./star_test

uninstall:
	@rm -f cubic_test
	@rm -f watts_strogatz_test
	@rm -f star_test