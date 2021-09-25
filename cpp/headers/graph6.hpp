#ifndef GRAPH6_H
#define GRAPH6_H
#include <string>
#include <wavefront.hpp>
#include <ip.hpp>
#include <genetic_algorithm.hpp>
using namespace std;
using namespace wf;
using namespace ip;
using namespace ga;
/* read graph6 for wavefront */
vector<wf::node> read_graph6_wf(const string &line);
/* read graph6 for ip */
ip::graph read_graph6_ip(const string &line);
/* read graph6 for genetic_algorithm */
ga::Graph read_graph6_ga(const string &line);

#endif