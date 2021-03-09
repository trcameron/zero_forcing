#include <iostream>
#include <fstream>
#include <vector>
#include <graph6.h>
#include <ilcplex/ilocplex.h>
using namespace std;

int zero_forcing(vector<vector<int>> adj)
{
	IloEnv env;
	/*
	IloModel model(env);
	IloNumVarArray var(env);
	IloRangeArray con(env);
	*/
	
	return 0;
}

int main(int argc,char** argv)
{
	/*
	const string line = "GruZp_";
	*/
	vector<vector<int>> adj;
	string line;
	ifstream file;
	if (argc > 1)
	{
		file.open(argv[1]);
	}
	istream &f = (argc > 1 ? file : cin);
	while(getline(f,line))
	{
		adj = read_graph6(line);
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