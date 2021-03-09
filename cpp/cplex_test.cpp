#include <iostream>
#include <fstream>
#include <vector>
#include <graph6.h>
#include <ilcplex/ilocplex.h>
using namespace std;

int zero_forcing(vector<vector<int>> adj)
{
	/* graph order */
	int n = adj.size();
	/* maximal propogation time */
	int T = n-1;
	/* graph edges */
	vector<pair<int,int>> edges;
	pair<int,int> e;
	int m = 0;
	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			if(adj[i][j]==1)
			{
				e.first = i;
				e.second = j;
				edges.push_back(e);
				m += 1;
			}
		}
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
		for(int w=0; w<n; w++)
		{
			if(w!=edges[k].second && adj[edges[k].first][w]==1)
			{
				con.add(IloRange(env,-IloInfinity,T));
				con[count].setLinearCoef(var[n+w],1);
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
       env.error() << "Failed to optimize LP" << endl;
       throw(-1);
    }
    IloNumArray vals(env);
    env.out() << "Solution status = " << cplex.getStatus() << endl;
    env.out() << "Solution value  = " << cplex.getObjValue() << endl;
	return n;
}

int main(int argc,char** argv)
{
	/*
	const string line = "GruZp_";
	*/
	vector<vector<int>> adj;
	string line;
	ifstream file;
	if(argc > 1)
	{
		file.open(argv[1]);
	}
	istream &f = (argc > 1 ? file : cin);
	while(getline(f,line))
	{
		adj = read_graph6(line);
		int n = zero_forcing(adj);
		/*int n = adj.size();*/
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<n; j++)
			{
				cout << adj[i][j] <<  " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	file.close();
	return 0;
}