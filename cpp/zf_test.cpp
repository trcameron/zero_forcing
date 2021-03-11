#include <iostream>
#include <fstream>
#include <vector>
#include <zero_forcing.h>
#include <ilcplex/ilocplex.h>
using namespace std;

int main(int argc,char** argv)
{
	vector<vector<int>> adj;
	string line;
	ifstream file;
	if(argc > 1)
	{
		file.open(argv[1]);
	}
	istream &f = (argc > 1 ? file : cin);
	while(getline(f,line))
	{
		adj = read_graph6(line);
		int zf = zf_ip(adj);
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
	}
	file.close();
	return 0;
}