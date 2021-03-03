#include <iostream>
#include <vector>
#include <graph6.h>
using namespace std;

vector<vector<int>> read_graph6(const string line)
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
	vector<vector<int>> adj(n,vector<int>(n,0));
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