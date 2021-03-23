#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <vector>
#include <stack>
#include <set>
#include <zero_forcing.h>
using namespace std;

/* Graph Closure */
void closure(graph g,set<int,less<int>> &s)
{
	/* set iterator */
	set<int,less<int>>::iterator it;
	/* get order and egdes of g*/
	int order = g.get_order();
	vector<pair<int,int>> edges = g.get_edges();
	/* initialize colored and count vectors and active nodes stack*/
	vector<bool> colored(order,0);
	vector<int> count(order,0);
	stack<int> active;
	for(it=s.begin(); it!=s.end(); it++)
	{
		colored[*it] = true;
	}
	for(it=s.begin(); it!=s.end(); it++)
	{
		vector<int> nbhd = g.get_neighbors(*it);
		for(int j=0; j<nbhd.size(); j++)
		{
			count[*it] += colored[nbhd[j]];
		}
		if(count[*it]==(g.get_degree(*it)-1))
		{
			active.push(*it);
		}
	}
	/* while there are still vertices that can do forcing*/
	while(!active.empty())
	{
		int u = active.top(); active.pop();		// vertex that forces
		vector<int> nbhd = g.get_neighbors(u);	// neighborhood of forcing vertex
		for(int i=0; i<nbhd.size(); i++)
		{
			if(!colored[nbhd[i]])				// vertex that is forced
			{
				colored[nbhd[i]] = true;		// color vertex
				s.insert(nbhd[i]);				// add to set for closure
				vector<int> nbhd2 = g.get_neighbors(nbhd[i]);	// neighborhood of forced vertex
				// update the number of colored neighbors for each colored vertex in nbhd2
				for(int j=0; j<nbhd2.size(); j++)
				{
					if(colored[nbhd2[j]])
					{
						count[nbhd2[j]] += 1;
						if(count[nbhd2[j]]==(g.get_degree(nbhd2[j])-1))
						{
							active.push(nbhd2[j]);
						}
					}
				}
				// count the number of colored neighbors of forced vertex
				for(int j=0; j<nbhd2.size(); j++)
				{
					count[nbhd[i]] += colored[nbhd2[j]];
				}
				if(count[nbhd[i]]==(g.get_degree(nbhd[i])-1))
				{
					active.push(nbhd[i]);
				}
				// break out of for i loop
				break;
			}
		}			
	}
}

/* Wavefront */
int wavefront(graph g)
{
	// graph order
	int order = g.get_order();
	// closure pairs
	vector<pair<set<int,less<int>>,int>> cl_pairs = {{{},0}};
	// iterate over all possible cardinalities
	for(int i=1; i<=order; i++)
	{
		// iterate over all closure pairs
		for(int j=0; j<cl_pairs.size(); j++)
		{
			// initialize set s and value r
			set<int,less<int>> s = cl_pairs[j].first;
			int r = cl_pairs[j].second;
			// iterate over all vertices
			for(int k=0; k<order; k++)
			{
				// initialize r_new and s_new
				int r_new = r;
				set<int,less<int>> s_new(s);
				// update r_new if k is not in s_new, also add k to s_new
				if(s_new.find(k)==s_new.end())
				{
					r_new += 1;
					s_new.insert(k);
				}
				//update r_new with max((nbhd-s_new).size()-1,0), also add nbhd to s_new
				int count = -1;
				vector<int> nbhd = g.get_neighbors(k);
				for(int l=0; l<nbhd.size(); l++)
				{
					if(s_new.find(nbhd[l])==s_new.end())
					{
						count += 1;
						s_new.insert(nbhd[l]);
					}
				}
				r_new += max(count,0);
				// compute set closure
				closure(g,s_new);
				// check if (s_new,r_new) needs to be added to cl_pairs
				if(r_new<=i)
				{
					pair<set<int,less<int>>,int> e = make_pair(s_new,0);
					while(find(cl_pairs.begin(),cl_pairs.end(),e)==cl_pairs.end() && e.second<=i)
					{
						e.second += 1;
					}
					if(e.second==(i+1))
					{
						cl_pairs.push_back({s_new,r_new});
						// check to return r
						if(s_new.size()==order)
						{
							return r_new;
						}
					}
				}
			}
		}
	}
	return 0;
}

/* main function */
int main(int argc,char** argv)
{
	/*
	graph g(4,{{0,2},{0,3},{1,3}});
	
	g.print_nodes();
	g.print_edges();
	
	vector<int> nbhd = g.get_neighbors(2);
	cout << "nbhd of 2: ";
	for(int i=0; i<nbhd.size(); i++)
	{
		cout << nbhd[i] << " ";
	}
	cout << endl;
	cout << "deg of 2: " << g.get_degree(2) << endl;
	
	set<int,less<int>> s = {2};
	set<int,less<int>>::iterator it;
	closure(g,s);
	cout << "closure of {2}: ";
	for(it=s.begin(); it!=s.end(); it++)
	{
		cout << *it << " ";
	}
	cout << endl;
	
	int zf = wavefront(g);
	cout << "zf number: " << zf << endl;
	return 0;
	*/
	
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
		graph g = read_graph6(line);
		int zf = wavefront(g);
		count += 1;
		g.print_nodes();
		g.print_edges();
		cout << "zf number: " << zf << "\n" << endl;
	}
	file.close();
	auto stop = chrono::high_resolution_clock::now();
	auto duration = duration_cast<chrono::microseconds>(stop - start);
	cout << "Average Elapsed Time of Wavefront Solver: "
	         << (duration.count()/count)*1E-6 << " seconds" << endl;
	return 0;
}