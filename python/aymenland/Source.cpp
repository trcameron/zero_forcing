#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <utility>
#include <list>
#include <fstream>
#include <bitset>
#include <unordered_map>
using namespace std;
using namespace std::chrono;

constexpr auto M = 300;


struct Bitset {
    bitset<M> bit;
    int size;

    Bitset(){}

    Bitset(int order) { size = order; }

    Bitset(vector<int>& v, int order) {
        size = order;
        for (int i = 0; i < v.size(); i++)
            bit.set(v[i]);
    }

    vector<int> convert() {
        vector<int> output;
        for (int i = 0; i < size; i++)
            if (bit[i])
                output.push_back(i);
        return output;
    }

    bool operator==(const Bitset& rhs) const {
        return bit == rhs.bit && size == rhs.size;
    }
};

Bitset operator| (const Bitset& lhs, const Bitset& rhs) {
    Bitset output(lhs.size);
    output.bit = lhs.bit | rhs.bit;
    return output;
}


void sample(vector<int>& vertices, Bitset& newVector, int k) {
    random_shuffle(vertices.begin(), vertices.end());
    for (int i = k; i--; )
        newVector.bit.set(vertices[i]);
}

void sample(vector<int>& vertices, vector<int>& newVector, int k) {
    random_shuffle(vertices.begin(), vertices.end());
    for (int i = k; i--; )
        newVector.push_back(vertices[i]);
}


/*bool In(int vertex, const vector<int>& vertices) {
    // Returns true if vertex is not in vertices, false otherwise
    for (int i = vertices.size(); i--; )
        if (vertex == vertices[i])
            return true;

    return false;
}*/


class Graph
{
protected:
    int order;
    vector<pair<int, int>> edges;
    vector<bool> colors;

public:
    vector< Bitset> neighbors;
    int max_degree;
    int min_degree;

    Graph() {}

    Graph(int order0, const vector<pair<int, int>>& edges0)
    {
        order = order0;
        for (int i = 0; i < edges0.size(); i++)
        {
            edges.push_back(edges0[i]);
        }

        neighbors = {};
        max_degree = 0;
        min_degree = 9999999;
        for (int j = 0; j < order; j++) {
            colors.push_back("white");
            Bitset nbhd(order);
            for (int i = 0; i < edges.size(); i++) {
                if (edges[i].first == j) {
                    nbhd.bit.set(edges[i].second);
                }
                else if (edges[i].second == j) {
                    nbhd.bit.set(edges[i].first);
                }
            }

            if (nbhd.bit.count() > max_degree)
                max_degree = nbhd.bit.count();
            if (nbhd.bit.count() < min_degree)
                min_degree = nbhd.bit.count();

            neighbors.push_back(nbhd);
        }
    }
    
    int get_order() const
    {
        return order;
    }

    vector<pair<int, int>> get_edges() const
    {
        return edges;
    }

    vector<int> vertices() const {
        vector<int> ret(order, 0);
        for (int i = 0; i < order; i++)
            ret[i] = i;

        return ret;
    }

    vector<vector<int>> get_adj() const
    {
        vector<vector<int>> adj(order, vector<int>(order, 0));
        for (int k = 0; k < edges.size(); k++)
        {
            adj[edges[k].first][edges[k].second] = 1;
            adj[edges[k].second][edges[k].first] = 1;
        }
        return adj;
    }

    Bitset adj(int node) const{
        return neighbors[node];
    }

    int get_degree(int node) const {
        return neighbors[node].bit.count();
    }

    void setColor(int node, const bool& color) {
        colors[node] = color;
    }

    bool getColor(int node) const {
        return colors[node];
    }

    void setAllColor(const bool& color) {
        for (int i = 0; i < order; i++)
            colors[i] = color;
    }
};


// Wavefront test
/*
void closure(Graph& g, set<int, less<int>>& s)
{
    // set iterator
    set<int, less<int>>::iterator it;
    // get order and egdes of g
    int order = g.get_order();
    vector<pair<int, int>> edges = g.get_edges();
    // initialize colored and count vectors and active nodes stack
    vector<bool> colored(order, 0);
    vector<int> count(order, 0);
    stack<int> active;
    for (it = s.begin(); it != s.end(); it++)
    {
        colored[*it] = true;
    }
    for (it = s.begin(); it != s.end(); it++)
    {
        vector<int> nbhd = g.adj(*it);
        for (int j = 0; j < nbhd.size(); j++)
        {
            count[*it] += colored[nbhd[j]];
        }
        if (count[*it] == (g.get_degree(*it) - 1))
        {
            active.push(*it);
        }
    }
    // while there are still vertices that can do forcing
    while (!active.empty())
    {
        int u = active.top(); active.pop();		// vertex that forces
        vector<int> nbhd = g.adj(u);			// neighborhood of forcing vertex
        for (int i = 0; i < nbhd.size(); i++)
        {
            if (!colored[nbhd[i]])				// vertex that is forced
            {
                colored[nbhd[i]] = true;		// color vertex
                s.insert(nbhd[i]);				// add to set for closure
                vector<int> nbhd2 = g.adj(nbhd[i]);	// neighborhood of forced vertex
                // update the number of colored neighbors for each colored vertex in nbhd2
                for (int j = 0; j < nbhd2.size(); j++)
                {
                    if (colored[nbhd2[j]])
                    {
                        count[nbhd2[j]] += 1;
                        if (count[nbhd2[j]] == (g.get_degree(nbhd2[j]) - 1))
                        {
                            active.push(nbhd2[j]);
                        }
                    }
                }
                // count the number of colored neighbors of forced vertex
                for (int j = 0; j < nbhd2.size(); j++)
                {
                    count[nbhd[i]] += colored[nbhd2[j]];
                }
                if (count[nbhd[i]] == (g.get_degree(nbhd[i]) - 1))
                {
                    active.push(nbhd[i]);
                }
                // break out of for i loop
                break;
            }
        }
    }
}

int zf_wave(Graph& g) {
    // graph order
    int order = g.get_order();
    // closure pairs
    vector<pair<set<int, less<int>>, int>> cl_pairs = { {{},0} };
    // iterate over all possible cardinalities
    for (int i = 1; i <= order; i++)
    {
        // iterate over all closure pairs
        for (int j = 0; j < cl_pairs.size(); j++)
        {
            // initialize set s and value r
            set<int, less<int>> s = cl_pairs[j].first;
            int r = cl_pairs[j].second;
            // iterate over all vertices
            for (int k = 0; k < order; k++)
            {
                // initialize r_new and s_new
                int r_new = r;
                set<int, less<int>> s_new(s);
                // update r_new if k is not in s_new, also add k to s_new
                if (s_new.find(k) == s_new.end())
                {
                    r_new += 1;
                    s_new.insert(k);
                }
                //update r_new with max((nbhd-s_new).size()-1,0), also add nbhd to s_new
                int count = -1;
                vector<int> nbhd = g.adj(k);
                for (int l = 0; l < nbhd.size(); l++)
                {
                    if (s_new.find(nbhd[l]) == s_new.end())
                    {
                        count += 1;
                        s_new.insert(nbhd[l]);
                    }
                }
                r_new += max(count, 0);
                // compute set closure
                closure(g, s_new);
                // check if (s_new,r_new) needs to be added to cl_pairs
                if (r_new <= i)
                {
                    pair<set<int, less<int>>, int> e = make_pair(s_new, 0);
                    while (find(cl_pairs.begin(), cl_pairs.end(), e) == cl_pairs.end() && e.second <= i)
                    {
                        e.second += 1;
                    }
                    if (e.second == (i + 1))
                    {
                        cl_pairs.push_back({ s_new,r_new });
                        // check to return r
                        if (s_new.size() == order)
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
*/

bool white_path(Graph& G, int s, int d) {
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


bool forcing_rule(Graph& G, int node0, const Bitset& b) {
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


bool psd_rule(Graph& G, int node0, const Bitset& b) {
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


struct returnPair {
    Bitset vertices;
    int propagation;
};
returnPair forcing_process(Graph& G, const Bitset& b, bool (*rule)(Graph&, int, const Bitset&) = &forcing_rule, int t = 1) {
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

// Genetic Algorithm Skeleton

int population_size = 9;
int num_of_elite_chrom = 1;
const int tournament_selection_size = 4;
double mutation_rate = 0.25;
int min_size, max_size;
unordered_map<bitset<M>, int> fitness_map;

class Chromosome {
public:
    Graph G;
    bool (*rule)(Graph&, int, const Bitset&);
    Bitset genes;
    int fitness;
    int t;
    bool throttling_num;

    Chromosome(Graph& G1, bool (*rule1)(Graph&, int, const Bitset&) = &forcing_rule, bool throttling_num1 = false, float p = 0.5) {
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

    Chromosome(Graph& G1, Bitset& genes1, bool (*rule1)(Graph&, int, const Bitset&) = &forcing_rule, bool throttling_num1 = false) {
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
            return -9999999;
        t = res.propagation;
        int result;
        if (throttling_num)
            result = 9999999 - t - genes.bit.count();
        else
            result = 9999999 - t - pow(genes.bit.count() + 2, 4);

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


class Population {
public:
    vector<Chromosome> chromosomes;

    Population(int size, Graph& G1, bool (*rule1)(Graph&, int, const Bitset&) = &forcing_rule, bool throttling_num = false) {
        for (int i = 0; i < size; i++) {
            chromosomes.push_back(Chromosome(G1, rule1, throttling_num, (i + 1) * 0.1));
        }
    }
};


class GeneticAlgorithm {
public:
    static Population evolve(Population& pop, Graph& G, bool (*rule)(Graph&, int, const Bitset&) = &forcing_rule, bool throttling_num = false) {
        Population pop1 = GeneticAlgorithm::crossover_population(pop, G, rule, throttling_num);
        GeneticAlgorithm::mutate_population(pop1, G);
        for (int i = 0; i < pop1.chromosomes.size(); i++)
            pop1.chromosomes[i].fitness = pop1.chromosomes[i].get_fitness();
        return pop1;
    }


    static Population crossover_population(const Population& pop, Graph& G, bool (*rule)(Graph&, int, const Bitset&) = &forcing_rule, bool throttling_num = false) {
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


    static Chromosome crossover_chromosomes(Bitset& genes1, Bitset& genes2, Graph& G, bool (*rule)(Graph&, int, const Bitset&) = &forcing_rule, bool throttling_num = false) {
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
        
        int top = -99999999;
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


void print_population(Population& pop, int gen_number) {
    cout << endl << "--------------------------------------------------" << endl;
    cout << "Generation #" << gen_number << "| Fittest chromosome fitness:" << pop.chromosomes[0].get_fitness() << endl;
    cout << "--------------------------------------------------" << endl;
    int i = 0;
    for (int i = 0; i < pop.chromosomes.size(); i++) {
        cout << "Chromosome  #" << i + 1 << " :" << pop.chromosomes[i] << "| Fitness: " << pop.chromosomes[i].get_fitness() << endl;
    }
}

// Applying Genetic Algorithm

struct returnTriplet {
    int zero_forcing_num;
    int propagation;
    vector<int> zero_forcing_set;
    int throttling_num;

    // int reached_target;
};
returnTriplet zero_forcing(Graph& G, bool (*rule)(Graph&, int, const Bitset&) = &forcing_rule, bool throttling_num = false) {
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

    return { (int)population.chromosomes[0].genes.bit.count(), population.chromosomes[0].t, population.chromosomes[0].genes.convert(), (int)population.chromosomes[0].genes.bit.count() + population.chromosomes[0].t }; //, 501 - reached_target};
}


int main() {
    
    // Data Gathering
    /*
    fstream file("file.txt");
    fstream out("output.txt");
    int n1(0), n2(0);
    int order = 0;

    while (file >> order){
        cout << order << endl;
        vector< pair<int, int>> edges;
        while (true) {
            file >> n1;
            if (n1 == -1)
                break;
            file >> n2;
            edges.push_back({ n1, n2 });
        }
        Graph G(order, edges);

        returnTriplet z = zero_forcing(G, forcing_rule);
        out << order << " " << G.get_edges().size() << " " << z.reached_target << endl;
    }*/
   
    // Test against Wavefront
    /*fstream file("graphs7.txt");
    int n1(0), n2(0);
    int order = 0;
    int fails = 0;

    while (file >> order) {
        vector< pair<int, int>> edges;
        while (true) {
            file >> n1;
            if (n1 == -1)
                break;
            file >> n2;
            edges.push_back({ n1, n2 });
        }
        if (edges.size() == 0)
            continue;

        Graph G(order, edges);

        auto start = high_resolution_clock::now();
        returnTriplet z = zero_forcing(G, forcing_rule);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        auto start1 = high_resolution_clock::now();
        int w = zf_wave(G);
        auto stop1 = high_resolution_clock::now();
        auto duration1 = duration_cast<microseconds>(stop1 - start1);

        cout << "Order: " << order << "\nn_Edges: " << G.get_edges().size() << "\nGA: " << z.zero_forcing_num << "\tTime: " << duration.count() / 1000000.0 << "\nWV: " << w << "\tTime: " << duration1.count() / 1000000.0 << endl << endl;
        if (z.zero_forcing_num != w)
            fails++;
    }
    cout << "Mismatch: " << fails << endl;*/
    
    Graph G(17, { {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 9}, {9, 10}, {10, 11}, {11, 12}, {12, 13}, {13, 14} });
    auto start = high_resolution_clock::now();
    returnTriplet z = zero_forcing(G, forcing_rule);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << duration.count() / 1000000.0 << endl;

    cout << "1. Zero Forcing number: " << z.zero_forcing_num << "  |  Propagation number: " << z.propagation << "  |  Zero Forcing Set: ";

    cout << "[";
    vector<int> lst = z.zero_forcing_set;

    for (int i = 0; i < lst.size(); i++) {
        if (i == lst.size() - 1)
            cout << lst[i];
        else
            cout << lst[i] << ", ";
    }
    cout << "]  |  Throttling Number: " << z.throttling_num << endl << /*"2. Wavefront: " << zf_wave(G) <<*/ endl;

    return 0;
}
