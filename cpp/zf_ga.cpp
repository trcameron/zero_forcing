#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <algorithm>
#include <random>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <zero_forcing.h>
using namespace std;
using namespace std::chrono;


void sample(const vector<int>& vertices, vector<int>& newVector, int k) {
    vector<int> v = vertices;
	random_device rd;
	mt19937 g(rd());
    std::shuffle(v.begin(), v.end(), g);
    for (int i = k; i--; )
        newVector.push_back(v[i]);
}


bool In(int vertex, const vector<int>& vertices) {
    /* Returns true if vertex is not in vertices, false otherwise */
    for (int i = vertices.size(); i--; )
        if (vertex == vertices[i])
            return true;

    return false;
}


class Graph {
    public:
    vector<string> colors;
    vector< vector<int>> adjacencies;
    vector<int> vertices;
    vector< vector<int>> edges;

    Graph(const vector< vector<int>>& edges1, bool adj_matrix = false) {
        if (!adj_matrix)
            edges = edges1;
        else
            for (int i = 0; i < edges1.size(); ++i)
                for (int j = i; j < edges1[0].size(); ++j)
                    if (edges1[i][j] != 0)
                        edges.push_back({ i, j });

        for (int i = 0; i < edges.size(); i++) {
            if (!In(edges[i][0], vertices)) {
                vertices.push_back(edges[i][0]);
                colors.push_back("white");
            }
            if (!In(edges[i][1], vertices)) {
                vertices.push_back(edges[i][1]);
                colors.push_back("white");
            }
        }

        for (int i = 0; i < vertices.size(); i++) {
            vector<int> lst;
            for (int j = 0; j < edges.size(); j++) {
                if (vertices[i] == edges[j][0]) {
                    if (!In(edges[j][1], lst))
                        lst.push_back(edges[j][1]);
                }
                else if (vertices[i] == edges[j][1])
                    if (!In(edges[j][0], lst))
                        lst.push_back(edges[j][0]);
            }

            adjacencies.push_back(lst);
        }
    }

    Graph(){}

    vector<int> adj(int vertex) {
		cout << "testing: " << vertices.size() << endl;
        for (int i = 0; i < vertices.size(); i++) {
            if (vertices[i] == vertex)
                return adjacencies[i];
        }
    }

    string getColor(int vertex) {
        for (int i = 0; i < vertices.size(); i++) {
            if (vertices[i] == vertex)
                return colors[i];
        }
    }

    void setColor(int vertex, const string& color) {
        for (int i = 0; i < vertices.size(); i++) {
            if (vertices[i] == vertex) {
                colors[i] = color;
                return;
            }
        }
    }

    void setAllColor(const string& color) {
        for (int i = 0; i < vertices.size(); i++)
            colors[i] = color;
    }
};


bool forcing_rule(Graph& G, int node0, const vector<int>& b) {
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


bool psd_rule(Graph& G, int node0, const vector<int>& b) {
    /* Returns True if, given a graph and a set of blue vertices b, nodes0 will be forced blue in the next iteration
    using the PSD forcing rule.
    Returns False otherwise. */

    if (G.getColor(node0) == "blue")
        return true;

    /* For every blue node node1 that is adjancent to node0, if node0 is not adjancent to any other white vertex adjacent to node1, then node0 will be forced in the next iteration, therefore return True. If this doesn't apply for any blue node, then node0 won't be forced blue, therefore return False. */
    vector<int> adj_node0 = G.adj(node0);
    for (int i = 0; i < b.size(); i++)
        if (In(b[i], adj_node0)) {
            bool works = true;
            vector<int> lst = G.adj(b[i]);
            for (int j = 0; j < lst.size(); j++) {
                if (In(lst[j], adj_node0) && G.getColor(lst[j]) == "white" && lst[j] != node0) {
                    works = false;
                    break;
                }
            }

            if (works)
                return true;
        }

    return false;
}


struct returnPair {
    vector<int> vertices;
    int propagation;
};
returnPair forcing_process(Graph& G, const vector<int>& b, bool (*rule)(Graph&, int, const vector<int>&) = &forcing_rule, int t = 1) {
    /* Applies the forcing rule given to the given graph with the given initial set of blue vertices b. Once the forcing process is done, returns the final set of blue vertices. */
    vector<int> b1 = b;
    while (true) {
        string blue = "blue";
        for (int i = 0; i < b1.size(); i++)
            G.setColor(b1[i], blue);


        if (b1.size() == G.vertices.size())
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

        if (b1 == new_b || new_b.size() == G.vertices.size())
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

class Chromosome {
public:
    Graph G;
    bool (*rule)(Graph&, int, const vector<int>&);
    vector<int> genes;
    int fitness;
    int t;
    bool throttling_num;

    Chromosome(Graph& G1, bool (*rule1)(Graph&, int, const vector<int>&) = &forcing_rule, bool throttling_num1 = false) {
        G = G1;
        rule = rule1;
        int k = rand() % G.vertices.size() + 1;
        sample(G.vertices, genes, k);
        t = 0;
        throttling_num = throttling_num1;
    }

    int get_fitness() {
        string white = "white";
        G.setAllColor(white);

        returnPair res = forcing_process(G, genes, rule);
        int forcing = res.vertices.size();
        if (forcing != G.vertices.size())
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


class Population {
public:
    vector<Chromosome> chromosomes;

    Population(int size, Graph& G1, bool (*rule1)(Graph&, int, const vector<int>&) = &forcing_rule, bool throttling_num = false) {
        for (int i = 0; i < size; i++)
            chromosomes.push_back(Chromosome(G1, rule1, throttling_num));
    }
};


class GeneticAlgorithm {
public:
    static Population evolve(Population& pop, Graph& G, bool (*rule)(Graph&, int, const vector<int>&) = &forcing_rule, bool throttling_num = false) {
        Population pop1 = GeneticAlgorithm::crossover_population(pop, G, rule, throttling_num);
        Population pop2 = GeneticAlgorithm::mutate_population(pop1, G);
        for (int i = 0; i < pop2.chromosomes.size(); i++)
            pop2.chromosomes[i].fitness = pop2.chromosomes[i].get_fitness();
        return pop2;
    }


    static Population crossover_population(const Population& pop, Graph& G, bool (*rule)(Graph&, int, const vector<int>&) = &forcing_rule, bool throttling_num = false) {
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


    static Population mutate_population(Population& pop, const Graph& G) {
        for (int i = num_of_elite_chrom; i < population_size; i++)
            GeneticAlgorithm::mutate_chromosome(pop.chromosomes[i], G);
        return pop;
    }


    static Chromosome crossover_chromosomes(const Chromosome& chromosome1, const Chromosome& chromosome2, Graph& G, bool (*rule)(Graph&, int, const vector<int>&) = &forcing_rule, bool throttling_num = false) {
        Chromosome crossover_chrom(G, rule, throttling_num);
        crossover_chrom.genes = {};
        int k = rand() % (1 + max(chromosome1.genes.size(), chromosome2.genes.size()) - min(chromosome1.genes.size(), chromosome2.genes.size())) + min(chromosome1.genes.size(), chromosome2.genes.size());

        vector<int> combined_genes = chromosome1.genes;
        combined_genes.insert(combined_genes.end(), chromosome2.genes.begin(), chromosome2.genes.end());

        sample(combined_genes, crossover_chrom.genes, k);
        sort(crossover_chrom.genes.begin(), crossover_chrom.genes.end());
        crossover_chrom.genes.erase(unique(crossover_chrom.genes.begin(), crossover_chrom.genes.end()), crossover_chrom.genes.end());

        return crossover_chrom;
    }


    static void mutate_chromosome(Chromosome& chromosome, const Graph& G) {
        if ((rand() % 1000) / 1000.0 < mutation_rate) {
            vector<int> lst = G.vertices;
            //lst.insert(lst.end(), chromosome.genes.begin(), chromosome.genes.end());   add this line to make the chromosome's genes have twice the chance of being picked
            int k = rand() % G.vertices.size() + 1;
            chromosome.genes = {};
            sample(lst, chromosome.genes, k);
            //sort(chromosome.genes.begin(), chromosome.genes.end());                   IF YOU UNCOMMENT THE lst.insert LINE THEN UNCOMMENT THESE TWO ASWELL
            //chromosome.genes.erase(unique(chromosome.genes.begin(), chromosome.genes.end()), chromosome.genes.end());
        }
    }


    static vector<Chromosome> select_tournament_population(const Population& pop, bool (*rule)(Graph&, int, const vector<int>&) = &forcing_rule) {
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


void print_population(Population& pop, int gen_number) {
    cout << endl << "--------------------------------------------------" << endl;
    cout << "Generation #" << gen_number << "| Fittest chromosome fitness:" << pop.chromosomes[0].get_fitness() << endl;
    cout << "--------------------------------------------------" << endl;
    int i = 0;
    for (int i = 0; i < pop.chromosomes.size(); i++) {
        cout << "Chromosome  #" << i + 1 << " :" << pop.chromosomes[i] << "| Fitness: " << pop.chromosomes[i].get_fitness() << endl;
    }
}

// Applied Genetic Algorithm

struct returnTriplet {
    int zero_forcing_num;
    int propagation;
    vector<int> zero_forcing_set;
    int throttling_num;
};
returnTriplet zero_forcing(Graph& G, bool (*rule)(Graph&, int, const vector<int>&) = &forcing_rule, bool throttling_num = false) {
    /* Applies a Genetic Algorithm to the given graph in order to find its minimal zero-forcing set, and therefore zero-forcing number. Returns the zero-forcing number, propagation number, and zero-forcing set */

    double target_gen = max(37.708025691857216 * G.vertices.size() + 0.6188752011422203 * G.edges.size() - 225.01828721571377, 30.0);
    srand(time(NULL));

    Population population(population_size, G, rule, throttling_num);
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
        population = GeneticAlgorithm::evolve(population, G, rule, throttling_num);
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


int main(int argc,char** argv)
{
	int count = 0;
	vector<vector<int>> adj;
	string line;
	ifstream file;
	if(argc > 1)
	{
		file.open(argv[1]);
	}
	istream &f = (argc > 1 ? file : cin);
	auto start = chrono::high_resolution_clock::now();
	while(getline(f,line))
	{
		adj = read_graph6(line);
		int n = adj.size();
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<n; j++)
			{
				cout << adj[i][j] <<  " ";
			}
			cout << endl;
		}
		cout << endl;
		
		Graph G(adj,true);
		cout << "created graph" << endl;
		
		vector<int> v = G.adj(0);
		cout << "graph adjacency" << endl;
		
		for(int i=0; i<v.size(); i++)
		{
			cout << v[i] << endl;
		}
		
		returnTriplet zf = zero_forcing(G,psd_rule);
		count += 1;
	}
	file.close();
	auto stop = chrono::high_resolution_clock::now();
	auto duration = duration_cast<chrono::microseconds>(stop - start);
	cout << "Average Elapsed Time of IP Solver: "
	         << (duration.count()/count)*1E-6 << " seconds" << endl;
	return 0;
}
