#include <iostream>
#include <vector>
#include <zero_forcing.h>
#include <ilcplex/ilocplex.h>
using namespace std;

graph read_graph6(const string line)
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
	// build graph edges
	vector<pair<int,int>> edges;
	int k = 0;
	for(int j=1; j<n; j++)
	{
		for(int i=0; i<j; i++)
		{
			if(bits[k])
			{
				edges.push_back(make_pair(i,j));
			}
			k++;
		}
	}
	
	// return graph
	return graph(n,edges);
}

int zf_ip(graph g)
{
	/* graph order */
	int n = g.get_order();
	/* maximal propogation time */
	int T = n-1;
	/* directional edges */
	vector<pair<int,int>> edges;
	vector<pair<int,int>> ge = g.get_edges();
	int m = 0;
	for(int k=0; k<ge.size(); k++)
	{
		edges.push_back(ge[k]);
		edges.push_back(make_pair(ge[k].second,ge[k].first));
		m += 2;
	}
	/* Ilo model, variables, constraints, and objective function*/
	IloEnv env;
	IloModel model(env);
	IloNumVarArray var(env);
	IloRangeArray con(env);
	IloObjective obj  = IloMinimize(env);
	/* set Ilo variables and objective function */
	int count = 0;
	for(int i=0; i<n; i++)
	{
		var.add(IloNumVar(env,0,1,IloNumVar::Int));
		obj.setLinearCoef(var[count],1);
		count += 1;
	}
	for(int i=0; i<n; i++)
	{
		var.add(IloNumVar(env,0,T,IloNumVar::Int));
		obj.setLinearCoef(var[count],0);
		count += 1;
	}
	for(int i=0; i<m; i++)
	{
		var.add(IloNumVar(env,0,1,IloNumVar::Int));
		obj.setLinearCoef(var[count],0);
		count += 1;
	}
	/* constraint set 1:
		s_{j} + sum(y_{e}) = 1, where e = (i,j)
	*/
	count = 0;
	for(int j=0; j<n; j++)
	{
		con.add(IloRange(env,1,1));
		con[count].setLinearCoef(var[j],1);
		for(int k=0; k<m; k++)
		{
			if(edges[k].second==j)
			{
				con[count].setLinearCoef(var[2*n+k],1);
			}
		}
		count += 1;
	}
	/* constraint set 2:
		x_{i} - x_{j} + (T+1)y_{e} <= T, where e = (i,j)
	*/
	for(int k=0; k<m; k++)
	{
		con.add(IloRange(env,-IloInfinity,T));
		con[count].setLinearCoef(var[n+edges[k].first],1);
		con[count].setLinearCoef(var[n+edges[k].second],-1);
		con[count].setLinearCoef(var[2*n+k],T+1);
		count += 1;
	}
	/* constraint set 3:
		x_{w} - x_{j} + (T+1)y_{e} <= T, where e = (i,j), w!=j, w~i
	*/
	for(int k=0; k<m; k++)
	{
		vector<int> nbhd = g.get_neighbors(edges[k].first);
		vector<int>::iterator w;
		for(w=nbhd.begin(); w<nbhd.end(); w++)
		{
			if(*w != edges[k].second)
			{
				con.add(IloRange(env,-IloInfinity,T));
				con[count].setLinearCoef(var[n+*w],1);
				con[count].setLinearCoef(var[n+edges[k].second],-1);
				con[count].setLinearCoef(var[2*n+k],T+1);
				count += 1;
			}
		}
	}
	/* add objective function and constraints to model*/
    model.add(obj);
    model.add(con);
    // Optimize the problem and obtain solution.
	IloCplex cplex(model);
	cplex.setOut(env.getNullStream());
    if(!cplex.solve())
	{
       env.out() << "Failed to optimize LP" << endl;
    }
	/*
    env.out() << "Solution status = " << cplex.getStatus() << endl;
    env.out() << "Solution value  = " << cplex.getObjValue() << endl;
	*/
	// return
	return round(cplex.getObjValue());
}