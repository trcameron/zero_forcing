#include <genetic_algorithm.hpp>
using namespace std;
using namespace ga;
/* Bitset operator */
Bitset operator| (const Bitset& lhs, const Bitset& rhs) {
    Bitset output(lhs.size);
    output.bit = lhs.bit | rhs.bit;
    return output;
}
/* sample functions */
void ga::sample(vector<int>& vertices, Bitset& newVector, int k) {
    random_shuffle(vertices.begin(), vertices.end());
    for (int i = k; i--; )
        newVector.bit.set(vertices[i]);
}
void ga::sample(vector<int>& vertices, vector<int>& newVector, int k) {
    random_shuffle(vertices.begin(), vertices.end());
    for (int i = k; i--; )
        newVector.push_back(vertices[i]);
}
/* white_path */
bool ga::white_path(Graph& G, int s, int d) {
    if (s == d)
        return true;

    int n = G.get_order();

    bool* visited = new bool[n];
    for (int i = 0; i < n; i++)
        visited[i] = false;

    list<int> queue;
    visited[s] = true;
    queue.push_back(s);

    vector<int>::iterator i;

    while (!queue.empty())
    {
        s = queue.front();
        queue.pop_front();

        // Get all adjacent vertices of the dequeued vertex s
        // If a adjacent has not been visited, then mark it visited
        // and enqueue it
        Bitset lst = G.neighbors[s];
        for (int i = 0; i < lst.size; i++) {
            if (!lst.bit[i])
                continue;

            if (G.getColor(i))
                continue;
            // If this adjacent node is the destination node, then
            // return true
            if (i == d)
                return true;

            // Else, continue to do BFS
            if (!visited[i])
            {
                visited[i] = true;
                queue.push_back(i);
            }
        }
    }

    // If BFS is complete without visiting d
    return false;
}
/* std_rule */
bool ga::std_rule(Graph& G, int node0, const Bitset& b) {
    /* Returns True if, given a graph and a set of blue vertices b, nodes0 will be forced blue in the next iteration. Returns False otherwise. */

    if (G.getColor(node0))
        return true;

    /* For every blue node node1 that is adjancent to node0, if node0 is the only white node adjacent to node1, then the forcing rule applies, and node0 will be forced blue in the next iteration, therefore return True. If this doesn't apply for any blue node then node0 won't be forced blue, therefore return False. */

    for (int i = 0; i < b.size; i++) {
        if (!b.bit[i])
            continue;
        Bitset node_neighbors = G.neighbors[node0];
        if (node_neighbors.bit[i]) {
            bool works = true;
            Bitset neighbors = G.neighbors[i];
            for (int j = 0; j < neighbors.size; j++) {
                if (!neighbors.bit[j])
                    continue;
                if (!G.getColor(j) && j != node0) {
                    works = false;
                    break;
                }
            }

            if (works)
                return true;
        }
    }

    return false;
}
/* psd_rule */
bool ga::psd_rule(Graph& G, int node0, const Bitset& b) {
    /* Returns True if, given a graph and a set of blue vertices b, nodes0 will be forced blue in the next iteration
    using the PSD forcing rule.
    Returns False otherwise. */

    if (G.getColor(node0))
        return true;

    /* For every blue node node1 that is adjancent to node0, if node0 is not adjancent to any other white vertex adjacent to node1, then node0 will be forced in the next iteration, therefore return True. If this doesn't apply for any blue node, then node0 won't be forced blue, therefore return False. */
    Bitset adj_node0 = G.neighbors[node0];
    for (int i = 0; i < b.size; i++) {
        if (!b.bit[i])
            continue;
        if (adj_node0.bit[i]) {
            bool works = true;
            Bitset lst = G.neighbors[i];
            for (int j = 0; j < lst.size; j++) {
                if (!lst.bit[j])
                    continue;
                if (!G.getColor(j) && j != node0 && white_path(G, node0, j)) {
                    works = false;
                    break;
                }
            }

            if (works)
                return true;
        }
    }

    return false;
}
/* forcing_process */
returnPair ga::forcing_process(Graph& G, const Bitset& b, bool (*rule)(Graph&, int, const Bitset&), int t) {
    /* Applies the forcing rule given to the given graph with the given initial set of blue vertices b. Once the forcing process is done, returns the final set of blue vertices. */
    Bitset b1 = b;
    while (true) {
        for (int i = 0; i < b1.size; i++) {
            if (!b1.bit[i])
                continue;
            G.setColor(i, true);
        }

        if (b1.bit.count() == G.get_order())
            return { b1, t - 1 };

        Bitset new_b = b1;
        bitset<M> explored(0);
        for (int i = 0; i < b1.size; i++) {
            if (!b1.bit[i])
                continue;
            explored.set(i);
        }

        for (int i = 0; i < b1.size; i++) {
            if (!b1.bit[i])
                continue;
            Bitset neighbors = G.neighbors[i];
            for (int j = 0; j < neighbors.size; j++) {
                if (!neighbors.bit[j])
                    continue;
                if (!explored[j]) {
                    explored.set(j);
                    if (rule(G, j, b1))
                        new_b.bit.set(j);
                }
            }
        }

        if (b1 == new_b || new_b.bit.count() == G.get_order())
            return { new_b, t };

        b1 = new_b;
        ++t;
    }
}
/* heuristic */
returnPair ga::heuristic(Graph& G, bool (*rule)(Graph&, int, const Bitset&)){
	Bitset blue(G.get_order()), new_blue(G.get_order());
	returnPair closure, new_closure;
	
	blue.bit.set(0);
	closure = forcing_process(G,blue,rule);

	while(closure.vertices.size < G.get_order()){
		for(int i=0; i<G.get_order(); i++){									// loop over all vertices
			if(!blue.bit[i]){												// if this vertex is not currently set
				new_blue.bit.set(i);
				new_closure = forcing_process(G,new_blue,rule);
			
				if(new_closure.vertices.size > closure.vertices.size){		// if adding this vertex results in a larger closure
					blue = new_blue;												// update blue to new_blue and closure to new_closure
					closure = new_closure;
				}
			}
		}
		new_blue = blue;													// reset new_blue to blue for next round of iterations
	}
	return {blue,closure.propagation};										// return blue and propagation time?
}
// Genetic Algorithm Skeleton
const int population_size = 9;
const int num_of_elite_chrom = 1;
const int tournament_selection_size = 4;
const double mutation_rate = 0.25;
int min_size, max_size;
unordered_map<bitset<M>, int> fitness_map;
/* Chromosome class */
class Chromosome {
public:
    Graph G;
    bool (*rule)(Graph&, int, const Bitset&);
    Bitset genes;
    int fitness;
    int t;
    bool throttling_num;

    Chromosome(Graph& G1, bool (*rule1)(Graph&, int, const Bitset&) = &std_rule, bool throttling_num1 = false, float p = 0.5) {
        G = G1;
        rule = rule1;

        genes.size = G.get_order();
        vector<int> v = G.vertices();
         for (int i = 0; i < G1.get_order(); i++) {
            if ((float)rand() / RAND_MAX < p)
                genes.bit.set(i);
         }
        
        t = 0;
        throttling_num = throttling_num1;
    }

    Chromosome(Graph& G1, Bitset& genes1, bool (*rule1)(Graph&, int, const Bitset&) = &std_rule, bool throttling_num1 = false) {
        G = G1;
        rule = rule1;
        genes = genes1;
        t = 0;
        throttling_num = throttling_num1;
    }

    int get_fitness() {
        if (fitness_map.find(genes.bit) != fitness_map.end())
            return fitness_map[genes.bit];

        G.setAllColor(false);

        returnPair res = forcing_process(G, genes, rule);
        int forcing = res.vertices.bit.count();
        if (forcing != G.get_order())
            return INT_MIN;
        t = res.propagation;
        int result;
        if (throttling_num)
            result = INT_MAX - t - genes.bit.count();
        else
            // result = 9999999 - t - pow(genes.bit.count() + 2, 4);
			result = INT_MAX - t - pow(genes.bit.count()+2,4);

        G.setAllColor(false);

        fitness_map[genes.bit] = result;
        return result;
    }

    friend ostream& operator<<(ostream& os, const Chromosome& chr) {
        cout << "[";
        for (int i = 0; i < chr.genes.size; i++) {
            if (!chr.genes.bit[i])
                continue;
            cout << i << ", ";
        }
        cout << "]";
        return os;
    }

};
/* Population class */
class Population {
public:
    vector<Chromosome> chromosomes;

    Population(int size, Graph& G1, bool (*rule1)(Graph&, int, const Bitset&) = &std_rule, bool throttling_num = false) {
        for (int i = 0; i < size; i++) {
            chromosomes.push_back(Chromosome(G1, rule1, throttling_num, (i + 1) * 0.1));
        }
    }
};
/* GeneticAlgorithm class */
class GeneticAlgorithm {
public:
    static Population evolve(Population& pop, Graph& G, bool (*rule)(Graph&, int, const Bitset&) = &std_rule, bool throttling_num = false) {
        Population pop1 = GeneticAlgorithm::crossover_population(pop, G, rule, throttling_num);
        GeneticAlgorithm::mutate_population(pop1, G);
        for (int i = 0; i < pop1.chromosomes.size(); i++)
            pop1.chromosomes[i].fitness = pop1.chromosomes[i].get_fitness();
        return pop1;
    }


    static Population crossover_population(const Population& pop, Graph& G, bool (*rule)(Graph&, int, const Bitset&) = &std_rule, bool throttling_num = false) {
        Population crossover_pop(0, G, rule, throttling_num);
        for (int i = 0; i < num_of_elite_chrom; i++)
            crossover_pop.chromosomes.push_back(pop.chromosomes[i]);
        int i = num_of_elite_chrom;

        while (i < population_size) {
            Bitset genes1 = GeneticAlgorithm::select_tournament_population(pop), genes2 = GeneticAlgorithm::select_tournament_population(pop);
            crossover_pop.chromosomes.push_back(GeneticAlgorithm::crossover_chromosomes(genes1, genes2, G, rule, throttling_num));
            i += 1;
        }
        return crossover_pop;
    }


    static void mutate_population(Population& pop, const Graph& G) {
        for (int i = num_of_elite_chrom; i < population_size; i++)
            GeneticAlgorithm::mutate_chromosome(pop.chromosomes[i], G);
    }


    static Chromosome crossover_chromosomes(Bitset& genes1, Bitset& genes2, Graph& G, bool (*rule)(Graph&, int, const Bitset&) = &std_rule, bool throttling_num = false) {
        vector<int> genes = {};
        int n1 = genes1.bit.count(), n2 = genes2.bit.count();
        int k = rand() % (1 + max(n1, n2) - min(n1, n2)) + min(n1, n2);
        
        vector<int> genes1_vec = genes1.convert();
        vector<int> genes2_vec = genes2.convert();
        vector<int> combined_genes = genes1_vec;
        combined_genes.insert(combined_genes.end(), genes2_vec.begin(), genes2_vec.end());
        sample(combined_genes, genes, k);
        sort(genes.begin(), genes.end());
        genes.erase(unique(genes.begin(), genes.end()), genes.end());

        Bitset genes0(genes, G.get_order());

        Chromosome crossover_chrom(G, genes0, rule, throttling_num);
        return crossover_chrom;
    }


    static void mutate_chromosome(Chromosome& chromosome, const Graph& G) {
        if ((rand() % 1000) / 1000.0 < mutation_rate) {
            vector<int> lst = G.vertices();
            //lst.insert(lst.end(), chromosome.genes.begin(), chromosome.genes.end());   add this line to make the chromosome's genes have twice the chance of being picked
            int k = rand() % (max_size - min_size + 1) + min_size;
            chromosome.genes = Bitset(G.get_order());
            sample(lst, chromosome.genes, k);
            //sort(chromosome.genes.begin(), chromosome.genes.end());                   IF YOU UNCOMMENT THE lst.insert LINE THEN UNCOMMENT THESE TWO ASWELL
            //chromosome.genes.erase(unique(chromosome.genes.begin(), chromosome.genes.end()), chromosome.genes.end());
        }
    }


    static Bitset select_tournament_population(const Population& pop) {
        int tournament_pop[tournament_selection_size];
        int i = 0;
        while (i < tournament_selection_size) {
            tournament_pop[i] = rand() % population_size;
            i += 1;
        }
        
        int top = INT_MIN;
        int top_chrom = 0;
        for (int i = 0; i < tournament_selection_size; i++) {
            if (pop.chromosomes[tournament_pop[i]].fitness > top) {
                top = pop.chromosomes[tournament_pop[i]].fitness;
                top_chrom = tournament_pop[i];
            }
        }

        return pop.chromosomes[top_chrom].genes;
    }
};
/* zero_forcing */
returnQuad ga::zero_forcing(Graph& G, bool (*rule)(Graph&, int, const Bitset&), bool throttling_num) {
    /* Applies a Genetic Algorithm to the given graph in order to find its minimal zero-forcing set, and therefore zero-forcing number. Returns the zero-forcing number, propagation number, and zero-forcing set */
    int n = G.get_order();
    int e = G.get_edges().size();
    
    if (e == 0)
        return {n , 0, G.vertices(), n};

    max_size = n;
    min_size = e / n;
    double target_gen = max(38.0 * n - 250, 30.0); // 500;
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

    /*int reached_target = 1;
    vector<int> oldgenes = population.chromosomes[0].genes;*/

    while (generation_number <= target_gen) {
        population = GeneticAlgorithm::evolve(population, G, rule, throttling_num);
        sort(population.chromosomes.begin(), population.chromosomes.end(),
            [](Chromosome& a, Chromosome& b) -> bool
            {
                return a.fitness > b.fitness;
            });
        //print_population(population, generation_number);
        generation_number++;

        /*reached_target++;
        if (population.chromosomes[0].genes != oldgenes)
            reached_target = 0;
        oldgenes = population.chromosomes[0].genes;*/
    }
    unordered_map<bitset<M>, int> new_map;
    fitness_map = new_map;

    return { (int)population.chromosomes[0].genes.bit.count(), population.chromosomes[0].t, population.chromosomes[0].genes.convert(), (int)population.chromosomes[0].genes.bit.count() + population.chromosomes[0].t }; //, 501 - reached_target};
}