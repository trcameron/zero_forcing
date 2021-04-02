#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <zero_forcing.h>
using namespace std;

int main(int argc,char** argv)
{
	int count = 0;
	string line;
	ifstream file;
	if(argc > 1)
	{
		file.open(argv[1]);
	}
	istream &f = (argc > 1 ? file : cin);
	auto start = chrono::high_resolution_clock::now();
	while(getline(f,line))
	{
		graph g = read_graph6(line);
		int zf = zf_ip(g);
		count += 1;
		/*
		g.print_nodes();
		g.print_edges();
		cout << "zf number: " << zf << "\n" << endl;
		*/
	}
	file.close();
	auto stop = chrono::high_resolution_clock::now();
	auto duration = duration_cast<chrono::microseconds>(stop - start);
	cout << "Average Elapsed Time of IP Solver: "
	         << (duration.count()/count)*1E-6 << " seconds" << endl;
	return 0;
}