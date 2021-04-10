#include <ilcplex/ilocplex.h>
#include <iostream>
#include <random>
#include <vector>
#include <set>
#include <zero_forcing.h>
using namespace std;

/* read graph6 */
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
/* Zero Forcing IP */
int zf_ip(graph& g)
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
		vector<int> nbhd = g.adj(edges[k].first);
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
	// return
	return round(cplex.getObjValue());
}
/* Graph Closure */
void closure(graph& g,set<int,less<int>> &s)
{
	/* set iterator */
	set<int,less<int>>::iterator it;
	/* get order and egdes of g*/
	int order = g.get_order();
	vector<pair<int,int>> edges = g.get_edges();
	/* initialize colored and count vectors and active nodes stack*/
	vector<bool> colored(order,0);
	vector<int> count(order,0);
	stack<int> active;
	for(it=s.begin(); it!=s.end(); it++)
	{
		colored[*it] = true;
	}
	for(it=s.begin(); it!=s.end(); it++)
	{
		vector<int> nbhd = g.adj(*it);
		for(int j=0; j<nbhd.size(); j++)
		{
			count[*it] += colored[nbhd[j]];
		}
		if(count[*it]==(g.get_degree(*it)-1))
		{
			active.push(*it);
		}
	}
	/* while there are still vertices that can do forcing*/
	while(!active.empty())
	{
		int u = active.top(); active.pop();		// vertex that forces
		vector<int> nbhd = g.adj(u);			// neighborhood of forcing vertex
		for(int i=0; i<nbhd.size(); i++)
		{
			if(!colored[nbhd[i]])				// vertex that is forced
			{
				colored[nbhd[i]] = true;		// color vertex
				s.insert(nbhd[i]);				// add to set for closure
				vector<int> nbhd2 = g.adj(nbhd[i]);	// neighborhood of forced vertex
				// update the number of colored neighbors for each colored vertex in nbhd2
				for(int j=0; j<nbhd2.size(); j++)
				{
					if(colored[nbhd2[j]])
					{
						count[nbhd2[j]] += 1;
						if(count[nbhd2[j]]==(g.get_degree(nbhd2[j])-1))
						{
							active.push(nbhd2[j]);
						}
					}
				}
				// count the number of colored neighbors of forced vertex
				for(int j=0; j<nbhd2.size(); j++)
				{
					count[nbhd[i]] += colored[nbhd2[j]];
				}
				if(count[nbhd[i]]==(g.get_degree(nbhd[i])-1))
				{
					active.push(nbhd[i]);
				}
				// break out of for i loop
				break;
			}
		}			
	}
}
/* Wavefront */
int zf_wave(graph& g)
{
	// graph order
	int order = g.get_order();
	// closure pairs
	vector<pair<set<int,less<int>>,int>> cl_pairs = {{{},0}};
	// iterate over all possible cardinalities
	for(int i=1; i<=order; i++)
	{
		// iterate over all closure pairs
		for(int j=0; j<cl_pairs.size(); j++)
		{
			// initialize set s and value r
			set<int,less<int>> s = cl_pairs[j].first;
			int r = cl_pairs[j].second;
			// iterate over all vertices
			for(int k=0; k<order; k++)
			{
				// initialize r_new and s_new
				int r_new = r;
				set<int,less<int>> s_new(s);
				// update r_new if k is not in s_new, also add k to s_new
				if(s_new.find(k)==s_new.end())
				{
					r_new += 1;
					s_new.insert(k);
				}
				//update r_new with max((nbhd-s_new).size()-1,0), also add nbhd to s_new
				int count = -1;
				vector<int> nbhd = g.adj(k);
				for(int l=0; l<nbhd.size(); l++)
				{
					if(s_new.find(nbhd[l])==s_new.end())
					{
						count += 1;
						s_new.insert(nbhd[l]);
					}
				}
				r_new += max(count,0);
				// compute set closure
				closure(g,s_new);
				// check if (s_new,r_new) needs to be added to cl_pairs
				if(r_new<=i)
				{
					pair<set<int,less<int>>,int> e = make_pair(s_new,0);
					while(find(cl_pairs.begin(),cl_pairs.end(),e)==cl_pairs.end() && e.second<=i)
					{
						e.second += 1;
					}
					if(e.second==(i+1))
					{
						cl_pairs.push_back({s_new,r_new});
						// check to return r
						if(s_new.size()==order)
						{
							return r_new;
						}
					}
				}
			}
		}
	}
	return 0;
}
/* Zero Forcing Heuristic */
set<int,less<int>> heuristic(graph& g)
{
	// initialize set iterator
	set<int,less<int>>::iterator it;
	// initialize colored set and closure
	set<int,less<int>> z, c;
	// while z is not a zero-forcing set
	while(c.size() != g.get_order())
	{
		// find vertex that maximizes closure
		int a;
		set<int,less<int>> c_new;
		for(int v=0; v<g.get_order(); v++)
		{
			if(c.find(v)==c.end())
			{
				set<int,less<int>> s(c);
				s.insert(v);
				closure(g,s);
				if(s.size() > c_new.size())
				{
					a = v;
					c_new = s;
				}
			}
		}
		// update closure and add vertex
		c = c_new;
		z.insert(a);
		// check if any vertices can be removed from coloring
		set<int,less<int>> z_new(z);
		for(it=z_new.begin(); it!=z_new.end(); it++)
		{
			z.erase(*it);
			closure(g,z);
			if(z.size() != g.get_order())
			{
				z.insert(*it);
			}
		}
	}
	// return
	return z;
}



// Sample
void sample(const vector<int>& vertices, vector<int>& newVector, int k) {
    vector<int> v = vertices;
	random_device rd;
	mt19937 g(rd());
	shuffle(v.begin(), v.end(), g);
    for (int i = k; i--; )
        newVector.push_back(v[i]);
}
// In
bool In(int vertex, const vector<int>& vertices) {
    /* Returns true if vertex is not in vertices, false otherwise */
    for (int i = vertices.size(); i--; )
        if (vertex == vertices[i])
            return true;

    return false;
}
// Standard Forcing Rule
bool forcing_rule(graph& G, int node0, const vector<int>& b) {
    /* Returns True if, given a graph and a set of blue vertices b, nodes0 will be forced blue in the next iteration. Returns False otherwise. */

    if (G.getColor(node0) == "blue")
        return true;

    /* For every blue node node1 that is adjancent to node0, if node0 is the only white node adjacent to node1, then the forcing rule applies, and node0 will be forced blue in the next iteration, therefore return True. If this doesn't apply for any blue node then node0 won't be forced blue, therefore return False. */
    vector<int> adj_node0 = G.adj(node0);
    for (int i = 0; i < b.size(); i++)
        if (In(b[i], adj_node0)) {
            bool works = true;
            vector<int> lst = G.adj(b[i]);
            for (int j = 0; j < lst.size(); j++)
                if (G.getColor(lst[j]) == "white" && lst[j] != node0) {
                    works = false;
                    break;
                }

            if (works)
                return true;
        }

    return false;
}
// returnPair Struct
struct returnPair {
    vector<int> vertices;
    int propagation;
};
// Forcing Process
returnPair forcing_process(graph& G, const vector<int>& b, bool (*rule)(graph&, int, const vector<int>&) = &forcing_rule, int t = 1) {
    /* Applies the forcing rule given to the given graph with the given initial set of blue vertices b. Once the forcing process is done, returns the final set of blue vertices. */
    vector<int> b1 = b;
    while (true) {
        string blue = "blue";
        for (int i = 0; i < b1.size(); i++)
            G.setColor(b1[i], blue);

        if (b1.size() == G.get_order())
            return { b1, t - 1 };

        vector<int> new_b = b1;
        vector<int> explored = b1;
        for (int i = 0; i < b1.size(); i++) {
            vector<int> lst = G.adj(b1[i]);
            for (int j = 0; j < lst.size(); j++) {
                if (!In(lst[j], explored)) {
                    explored.push_back(lst[j]);
                    if (rule(G, lst[j], b1))
                        new_b.push_back(lst[j]);
                }
            }
        }

        if (b1 == new_b || new_b.size() == G.get_order())
            return { new_b, t };

        b1 = new_b;
        ++t;
    }
}
// Genetic Algorithm Skeleton
int population_size = 8;
int num_of_elite_chrom = 1;
int tournament_selection_size = 4;
double mutation_rate = 0.25;
int min_size, max_size;
// Chromosome Class
class Chromosome {
public:
    graph G;
    bool (*rule)(graph&, int, const vector<int>&);
    vector<int> genes;
    int fitness;
    int t;
    bool throttling_num;

    Chromosome(graph& G1, bool (*rule1)(graph&, int, const vector<int>&) = &forcing_rule, bool throttling_num1 = false) {
        G = G1;
        rule = rule1;
        int k = rand() % (max_size - min_size + 1) + min_size;
        sample(G.vertices(), genes, k);
        t = 0;
        throttling_num = throttling_num1;
    }

    int get_fitness() {
        string white = "white";
        G.setAllColor(white);

        returnPair res = forcing_process(G, genes, rule);
        int forcing = res.vertices.size();
        if (forcing != G.get_order())
            return -9999999;
        t = res.propagation;
        int result;
        if (throttling_num)
            result = 9999999 - t - genes.size();
        else
            result = 9999999 - t - pow(genes.size() + 2, 4);

        G.setAllColor(white);

        return result;
    }

    friend ostream& operator<<(ostream& os, const Chromosome& chr) {
        cout << "[";
        for (int i = 0; i < chr.genes.size(); i++) {
            if (i == chr.genes.size() - 1)
                cout << chr.genes[i];
            else
                cout << chr.genes[i] << ", ";
        }
        cout << "]";
        return os;
    }

};
// Population Class
class Population {
public:
    vector<Chromosome> chromosomes;

    Population(int size, graph& G1, bool (*rule1)(graph&, int, const vector<int>&) = &forcing_rule, bool throttling_num = false) {
        for (int i = 0; i < size; i++)
            chromosomes.push_back(Chromosome(G1, rule1, throttling_num));
    }
};
// Genetic Algorithm Class
class GeneticAlgorithm {
public:
    static Population evolve(Population& pop, graph& G, bool (*rule)(graph&, int, const vector<int>&) = &forcing_rule, bool throttling_num = false) {
        Population pop1 = GeneticAlgorithm::crossover_population(pop, G, rule, throttling_num);
        Population pop2 = GeneticAlgorithm::mutate_population(pop1, G);
        for (int i = 0; i < pop2.chromosomes.size(); i++)
            pop2.chromosomes[i].fitness = pop2.chromosomes[i].get_fitness();
        return pop2;
    }


    static Population crossover_population(const Population& pop, graph& G, bool (*rule)(graph&, int, const vector<int>&) = &forcing_rule, bool throttling_num = false) {
        Population crossover_pop(0, G, rule, throttling_num);
        for (int i = 0; i < num_of_elite_chrom; i++)
            crossover_pop.chromosomes.push_back(pop.chromosomes[i]);
        int i = num_of_elite_chrom;

        while (i < population_size) {
            Chromosome chromosome1 = GeneticAlgorithm::select_tournament_population(pop, rule)[0];
            Chromosome chromosome2 = GeneticAlgorithm::select_tournament_population(pop, rule)[0];
            crossover_pop.chromosomes.push_back(GeneticAlgorithm::crossover_chromosomes(chromosome1, chromosome2, G, rule, throttling_num));
            i += 1;
        }
        return crossover_pop;
    }


    static Population mutate_population(Population& pop, const graph& G) {
        for (int i = num_of_elite_chrom; i < population_size; i++)
            GeneticAlgorithm::mutate_chromosome(pop.chromosomes[i], G);
        return pop;
    }


    static Chromosome crossover_chromosomes(const Chromosome& chromosome1, const Chromosome& chromosome2, graph& G, bool (*rule)(graph&, int, const vector<int>&) = &forcing_rule, bool throttling_num = false) {
        Chromosome crossover_chrom(G, rule, throttling_num);
        crossover_chrom.genes = {};
        int n1 = chromosome1.genes.size(), n2 = chromosome2.genes.size();
        int k = rand() % (1 + max(n1, n2) - min(n1, n2)) + min(n1, n2);

        vector<int> combined_genes = chromosome1.genes;
        combined_genes.insert(combined_genes.end(), chromosome2.genes.begin(), chromosome2.genes.end());

        sample(combined_genes, crossover_chrom.genes, k);
        sort(crossover_chrom.genes.begin(), crossover_chrom.genes.end());
        crossover_chrom.genes.erase(unique(crossover_chrom.genes.begin(), crossover_chrom.genes.end()), crossover_chrom.genes.end());

        return crossover_chrom;
    }


    static void mutate_chromosome(Chromosome& chromosome, const graph& G) {
        if ((rand() % 1000) / 1000.0 < mutation_rate) {
            vector<int> lst = G.vertices();
            //lst.insert(lst.end(), chromosome.genes.begin(), chromosome.genes.end());   add this line to make the chromosome's genes have twice the chance of being picked
            int n = lst.size();
            int k = rand() % (max_size - min_size + 1) + min_size;
            chromosome.genes = {};
            sample(lst, chromosome.genes, k);
            //sort(chromosome.genes.begin(), chromosome.genes.end());                   IF YOU UNCOMMENT THE lst.insert LINE THEN UNCOMMENT THESE TWO ASWELL
            //chromosome.genes.erase(unique(chromosome.genes.begin(), chromosome.genes.end()), chromosome.genes.end());
        }
    }


    static vector<Chromosome> select_tournament_population(const Population& pop, bool (*rule)(graph&, int, const vector<int>&) = &forcing_rule) {
        vector<Chromosome> tournament_pop;
        int i = 0;
        while (i < tournament_selection_size) {
            tournament_pop.push_back(pop.chromosomes[rand() % population_size]);
            i += 1;
        }
        sort(tournament_pop.begin(), tournament_pop.end(),
            [](const Chromosome& a, const Chromosome& b) -> bool
            {
                return a.fitness > b.fitness;
            });

        return tournament_pop;
    }
};
// Zero Forcing GA
returnTriplet zf_ga(graph& G) {
    /* Applies a Genetic Algorithm to the given graph in order to find its minimal zero-forcing set, and therefore zero-forcing number. Returns the zero-forcing number, propagation number, and zero-forcing set */
    int n = G.get_order();
    int e = G.get_edges().size();
    
    if (e == 0)
        return {n , 0, G.vertices(), n};

    max_size = n;
    min_size = e / n;
    double target_gen = max(37.708025691857216 * n + 0.6188752011422203 * e - 255.01828721571377, 30.0);
    srand(time(NULL));

    Population population(population_size, G, forcing_rule, false);
    for (int i = 0; i < population.chromosomes.size(); i++)
        population.chromosomes[i].fitness = population.chromosomes[i].get_fitness();
    sort(population.chromosomes.begin(), population.chromosomes.end(),
        [](Chromosome& a, Chromosome& b) -> bool
        {
            return a.fitness > b.fitness;
        });
    //print_population(population, 0);
    int generation_number = 1;

    while (generation_number <= target_gen) {
        population = GeneticAlgorithm::evolve(population, G, forcing_rule, false);
        sort(population.chromosomes.begin(), population.chromosomes.end(),
            [](Chromosome& a, Chromosome& b) -> bool
            {
                return a.fitness > b.fitness;
            });
        //print_population(population, generation_number);
        generation_number++;
    }

    return { (int)population.chromosomes[0].genes.size(), population.chromosomes[0].t, population.chromosomes[0].genes, (int)population.chromosomes[0].genes.size() + population.chromosomes[0].t };
}