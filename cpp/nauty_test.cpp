#include <chrono>
#include <iostream>
#include <graph6.hpp>
#include <wavefront.hpp>
#include <genetic_algorithm.hpp>
using namespace std;
using namespace wf;
using namespace ga;
/* main function */
int main(int argc,char** argv)
{
	// open file
	ofstream fout;
	fout.open("data_files/nauty_test.dat");
	fout << "wf_time, ga_time, correct, error" << endl;
	// set precision
	fout.precision(5);
	// timing variables
	std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	// testing variables
	int corr = 0, count = 0;
	double err = 0, wf_duration = 0, ga_duration = 0;
	// input stream variables
	string line;
	ifstream file;
	if(argc > 1)
	{
		file.open(argv[1]);
	}
	istream &f = (argc > 1 ? file : cin);
	while(getline(f,line))
	{
		// wavefront
		start = std::chrono::high_resolution_clock::now();
		std::vector<wf::node> nodes = read_graph6_wf(line);
		int wf_zf_num = wf::zero_forcing(nodes);
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		wf_duration += duration.count()*1E-6;
		// genetic_algorithm
		start = std::chrono::high_resolution_clock::now();
		ga::Graph ga_graph = read_graph6_ga(line);
		returnQuad z = ga::zero_forcing(ga_graph,std_rule,false);
		int ga_zf_num = z.zero_forcing_num;
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		ga_duration += duration.count()*1E-6;
		// correctness
		if(ga_zf_num == wf_zf_num)
		{
			corr += 1;
		}
		// error in ga result
		err += abs(ga_zf_num - wf_zf_num);
		// count
		count += 1;
	}
	// write to file
	fout << std::scientific << wf_duration/count;
	fout << ", " << std::scientific << ga_duration/count;
	fout << ", " << std::fixed << (double)corr/count;
	fout << ", " << std::scientific << err/count << endl;
	// close file
	fout.close();
	// return
	return 0;
}