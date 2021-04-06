#include <chrono>
#include <iostream>
#include <fstream>
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
	vector<int> time_ip, time_wave, time_ga;
	while(getline(f,line))
	{
		graph g = read_graph6(line);
		
		auto start = chrono::high_resolution_clock::now();
		//int zf1 = zf_ip(g);
		auto stop = chrono::high_resolution_clock::now();
		auto duration = duration_cast<chrono::microseconds>(stop - start);
		time_ip.push_back(duration.count());
		
		start = chrono::high_resolution_clock::now();
		int zf2 = zf_wave(g);
		stop = chrono::high_resolution_clock::now();
		duration = duration_cast<chrono::microseconds>(stop - start);
		time_wave.push_back(duration.count());
		
		start = chrono::high_resolution_clock::now();
		returnTriplet z = zf_ga(g);
		stop = chrono::high_resolution_clock::now();
		duration = duration_cast<chrono::microseconds>(stop - start);
		time_ga.push_back(duration.count());
	}
	
	file.close();
	
	double avgtime_ip=0, avgtime_wave=0, avgtime_ga=0;
	for(int i=0; i<time_ip.size(); i++)
	{
		avgtime_ip += time_ip[i]*1E-6;
		avgtime_wave += time_wave[i]*1E-6;
		avgtime_ga += time_ga[i]*1E-6;
	}
	avgtime_ip = avgtime_ip/time_ip.size();
	avgtime_wave = avgtime_wave/time_ip.size();
	avgtime_ga = avgtime_ga/time_ip.size();

	cout << scientific << setprecision(3) << avgtime_ip << ", " << avgtime_wave << ", " << avgtime_ga << endl;
	
	return 0;
}