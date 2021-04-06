#include <iostream>
#include <fstream>
#include <zero_forcing.h>
using namespace std;

int main(int argc,char** argv)
{
	string line;
	ifstream file;
	if(argc > 1)
	{
		file.open(argv[1]);
	}
	istream &f = (argc > 1 ? file : cin);
	int count =0, corr_wave = 0, corr_ga = 0;
	while(getline(f,line))
	{
		graph g = read_graph6(line);
		int zf1 = zf_ip(g);
		int zf2 = zf_wave(g);
		returnTriplet z = zf_ga(g);
		int zf3 = z.zero_forcing_num;
		if(zf2 == zf1)
		{
			corr_wave += 1;
		}
		if(zf3 == zf1)
		{
			corr_ga += 1;
		}
		count += 1;
	}
	file.close();
	
	double pcorr_wave = static_cast<double>(corr_wave)/static_cast<double>(count);
	double pcorr_ga = static_cast<double>(corr_ga)/static_cast<double>(count);
	
	cout << setprecision(3) << pcorr_wave << ", " << pcorr_ga << endl;
	
	return 0;
}