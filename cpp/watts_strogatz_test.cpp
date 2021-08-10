#include <chrono>
#include <iostream>
#include <edg.hpp>
#include <wavefront.hpp>
#include <genetic_algorithm.hpp>
using namespace std;
using namespace wf;
using namespace ga;
/* main function */
int main(int argc,char** argv)
{
	// open file 5.3
	ofstream fout;
	fout.open("data_files/watts_strogatz_5.3_test.dat");
	fout << "order, wf_time, ga_time, correct, error" << endl;
	// set precision
	fout.precision(5);
	// timing variables
	std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	for(int i=10; i<=80; i+=10)
	{
		int corr = 0;
		double err = 0, wf_duration = 0, ip_duration = 0, ga_duration = 0;
		for(int j=1; j<=5; j++)
		{
			// file_name for graph in edg form
			string file_name = "fast_test_graphs/watts_strogatz_graphs/WS_"+to_string(i)+"_5_0.3_"+to_string(j)+".edg";
			// wavefront
			start = std::chrono::high_resolution_clock::now();
			std::vector<wf::node> nodes = read_edg_wf(file_name);
			int wf_zf_num = wf::zero_forcing(nodes);
			stop = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
			wf_duration += duration.count()*1E-6;
			// genetic_algorithm
			start = std::chrono::high_resolution_clock::now();
			ga::Graph ga_graph = read_edg_ga(file_name);
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
		}
		// write to file
		fout << i << ", " << std::scientific << wf_duration/5;
		fout << ", " << std::scientific << ga_duration/5;
		fout << ", " << std::fixed << (double)corr/5;
		fout << ", " << std::scientific << err/5 << endl;
	}
	cout << endl;
	// close file
	fout.close();
	// open file 10.3
	fout.open("data_files/watts_strogatz_10.3_test.dat");
	fout << "order, wf_time, ga_time, correct, error" << endl;
	// set precision
	fout.precision(5);
	for(int i=20; i<=70; i+=10)
	{
		int corr = 0;
		double err = 0, wf_duration = 0, ip_duration = 0, ga_duration = 0;
		for(int j=1; j<=5; j++)
		{
			// file_name for graph in edg form
			string file_name = "fast_test_graphs/watts_strogatz_graphs/WS_"+to_string(i)+"_10_0.3_"+to_string(j)+".edg";
			// wavefront
			start = std::chrono::high_resolution_clock::now();
			std::vector<wf::node> nodes = read_edg_wf(file_name);
			int wf_zf_num = wf::zero_forcing(nodes);
			stop = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
			wf_duration += duration.count()*1E-6;
			// genetic_algorithm
			start = std::chrono::high_resolution_clock::now();
			ga::Graph ga_graph = read_edg_ga(file_name);
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
		}
		// write to file
		fout << i << ", " << std::scientific << wf_duration/5;
		fout << ", " << std::scientific << ga_duration/5;
		fout << ", " << std::fixed << (double)corr/5;
		fout << ", " << std::scientific << err/5 << endl;
	}
	// close file
	fout.close();
}