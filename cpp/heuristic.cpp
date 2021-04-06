#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <numeric>
#include <zero_forcing.h>
using namespace std;

/* main function */
int main(int argc,char** argv)
{	
	string line;
	ifstream file;
	if(argc > 1)
	{
		file.open(argv[1]);
	}
	istream &f = (argc > 1 ? file : cin);
	vector<double> ratio;
	while(getline(f,line))
	{
		graph g = read_graph6(line);
		int zf_num = zf_wave(g);
		set<int,less<int>> zf_set = heuristic(g);
		ratio.push_back(static_cast<double>(zf_num)/static_cast<double>(zf_set.size()));
	}
	cout << setprecision(3) << *min_element(ratio.begin(),ratio.end()) << ", " << accumulate(ratio.begin(),ratio.end(),0.0)/ratio.size() << ", " << *max_element(ratio.begin(),ratio.end()) << endl;
	file.close();
	return 0;
}