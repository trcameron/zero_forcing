#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <zero_forcing.h>
#include <ilcplex/ilocplex.h>
using namespace std;

int main(int argc,char** argv)
{
	int count = 0;
	vector<vector<int>> adj;
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
		adj = read_graph6(line);
		int zf = zf_ip(adj);
		count += 1;
		/*
		int n = adj.size();
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<n; j++)
			{
				cout << adj[i][j] <<  " ";
			}
			cout << endl;
		}
		cout << endl;
		*/
	}
	file.close();
	auto stop = chrono::high_resolution_clock::now();
	auto duration = duration_cast<chrono::microseconds>(stop - start);
	cout << "Average Elapsed Time of IP Solver: "
	         << (duration.count()/count)*1E-6 << " seconds" << endl;
	return 0;
}