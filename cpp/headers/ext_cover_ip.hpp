#ifndef EXT_COVER_IP_H
#define EXT_COVER_IP_H
#include <ip.hpp>
namespace ecip
{
	/* cut_callback */
	int __stdcall cut_callback(GRBmodel *model, void *cbdata, int where, void *usrdata);
	/* zero_forcing */
	int zero_forcing(ip::graph& our_graph);
}
#endif