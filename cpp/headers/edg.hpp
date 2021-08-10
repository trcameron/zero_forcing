#ifndef EDG_H
#define EDG_H
#include <string>
#include <wavefront.hpp>
#include <ip.hpp>
#include <genetic_algorithm.hpp>
using namespace std;
using namespace wf;
using namespace ip;
using namespace ga;
/* read edg for wavefront */
std::vector<wf::node> read_edg_wf(const std::string &file_name);
/* read edg for ip */
ip::graph read_edg_ip(const std::string &file_name);
/* read edg for genetic_algorithm */
ga::Graph read_edg_ga(const std::string &file_name);

#endif