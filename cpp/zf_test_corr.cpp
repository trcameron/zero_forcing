#include <iostream>
#include <fstream>
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
	while(getline(f,line))
	{
		graph g = read_graph6(line);
		int zf1 = zf_ip(g);
		int zf2 = zf_wave(g);
		returnTriplet z = zf_ga(g);
		int zf3 = z.zero_forcing_num;
		if(zf1 != zf2)
		{
			cout << "IP and Wavefront disagree." << endl;
			cout << "ZF_IP = " << zf1 << endl;
			cout << "ZF_Wave = " << zf2 << endl;
		}
		if(zf1 != zf3)
		{
			cout << "IP and GA disagree." << endl;
			cout << "ZF_IP = " << zf1 << endl;
			cout << "ZF_GA = " << zf3 << endl;	
		}
		count += 1;
	}
	file.close();
	return 0;
}