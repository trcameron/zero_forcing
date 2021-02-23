#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

std::vector<std::vector<int>> read_graph6(const string line)
{
	int m = line.size();		// length of line
	int data[m];				// data for graph
	// shift line characters by 63
	for(int i=0; i<m; i++)
	{
		data[i] = line[i] - 63;	
	}
	// size and edge data
	int n;
	int* edge;
	if(data[0]<=62)
	{
		n = data[0];
		edge = &data[1];
		m -= 1;
	}
	else if(data[1]<=62)
	{
		n = (data[1]<<12) + (data[2]<<6) + data[3];
		edge = &data[4];
		m -= 4;
	}
	else
	{
		n = (data[2]<<30) + (data[3]<<24) + (data[4]<<18) + (data[5]<<12) + (data[6]<<6) + data[7];
		edge = &data[8];
		m -= 8;
	}
	// bits
	bool bits[m*6];
	for(int i=0;i<m;i++)
	{
		for(int j=5; j>=0; j--)
		{
			bits[i*6+(5-j)] = (edge[i]>>j) & 1;
		}
	}
	// initialize adjacency matrix
	std::vector<std::vector<int>> adj(n,std::vector<int>(n,0));
	int k = 0;
	for(int j=1; j<n; j++)
	{
		for(int i=0; i<j; i++)
		{
			if(bits[k])
			{
				adj[i][j] = 1;
				adj[j][i] = 1;
			}
			k++;
		}
	}
	return adj;
}

int main(int argc,char** argv)
{
	/*
	const string line = "GruZp_";
	*/
	std::vector<std::vector<int>> adj;
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
				std::cout << adj[i][j] <<  " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	file.close();
	return 0;
}