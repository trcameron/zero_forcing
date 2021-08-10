#include <chrono>
#include <iostream>
#include <edg.hpp>
#include <ip.hpp>
#include <infection_ip.hpp>
#include <genetic_algorithm.hpp>
using namespace std;
using namespace ifip;
using namespace ga;
/* main function */
int main(int argc,char** argv)
{
	// open file
	ofstream fout;
	fout.open("data_files/star_test.dat");
	fout << "order, ifip_time, ga_time, error" << endl;
	// set precision
	fout.precision(5);
	// timing variables
	std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	double err, ifip_duration, ga_duration;
	for(int i=10; i<=100; i+=10)
	{
		// file_name for graph in edg form
		string file_name = "fast_test_graphs/star_graphs/Star"+to_string(i)+".edg";
		// infection_ip
		start = std::chrono::high_resolution_clock::now();
		ip::graph our_graph = read_edg_ip(file_name);
		int ifip_zf_num = ifip::zero_forcing(our_graph);
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		ifip_duration = duration.count()*1E-6;
		// genetic_algorithm
		start = std::chrono::high_resolution_clock::now();
		ga::Graph ga_graph = read_edg_ga(file_name);
		returnQuad z = ga::zero_forcing(ga_graph,std_rule,false);
		int ga_zf_num = z.zero_forcing_num;
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		ga_duration = duration.count()*1E-6;
		// error in ga result
		err = ga_zf_num - ifip_zf_num;
		// write to file
		fout << (i+1) << ", " << std::scientific << ifip_duration;
		fout << ", " << std::scientific << ga_duration;
		fout << ", " << std::fixed << err << endl;
	}
	// close file
	fout.close();
}