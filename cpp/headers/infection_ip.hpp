#ifndef INFECTION_IP_H
#define INFECTION_IP_H
#include <ip.hpp>
namespace ifip
{
	/* cut_callback */
	int __stdcall cut_callback(GRBmodel *model, void *cbdata, int where, void *usrdata);
	/* zero_forcing */
	int zero_forcing(ip::graph& our_graph);
}
#endif