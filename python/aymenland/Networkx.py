import networkx as nx
import random


class GraphMatrix:
    """ Takes in a networkx graph and a list of diagonal entries starting from the top left to the bottom right, and
    creates the matrix corresponding to that graph. """

    def __init__(self, graph: nx.Graph, diagonal_entries: list):
        self.matrix = []

        c = 0
        for i in list(graph.nodes()):
            self.matrix.append([])
            for j in list(graph.nodes()):
                if i == j:
                    self.matrix[-1].append(diagonal_entries[c])
                else:
                    self.matrix[-1].append(1 if j in list(graph.adj[i]) else 0)

            c += 1

    """ Shows the matrix of the graph the object was initialized with. If upper is given to be True, shows the 
    upper-triangular form of that matrix. """

    def show(self, upper=False):
        print("[")
        for line in (self.matrix if not upper else self.upper_triangle()):
            print("[", end="")
            for i in range(len(line)):
                if i == len(line) - 1:
                    print((" %.2f" if line[i] >= 0 else "%.2f") % line[i], end="")
                else:
                    print((" %.2f" if line[i] >= 0 else "%.2f") % line[i] + ", ", end="")
            print("]")
        print("]")

    """ Forms and returns the upper-triangular form of the matrix of the given graph. """

    def upper_triangle(self):
        n = len(self.matrix)

        matrix = []
        for line in self.matrix:
            matrix.append(line.copy())

        for c in range(n - 1):
            i = c
            val = -1
            for j in range(len(matrix[c:])):
                if abs(matrix[c:][j][c]) > val:
                    val = abs(matrix[c:][j][c])
                    i = j + c

            if i != c:
                matrix[c], matrix[i] = matrix[i], matrix[c]

            if matrix[c][c] == 0:
                continue

            for line in matrix[c + 1:]:
                m = line[c] / matrix[c][c]
                for k in range(n):
                    line[k] -= m * matrix[c][k]

        return matrix

    """ Returns the rank of the matrix of the given graph. """

    def rank(self):
        n = len(self.matrix)
        return n - self.upper_triangle().count([0] * n)


def forcing_rule(graph: nx.Graph, node0, b):
    """ Returns True if, given a graph and a set of blue vertices b, nodes0 will be forced blue in the next iteration.
    Returns False otherwise. """

    if graph.nodes[node0]["color"] == "blue":
        return True

    """ For every blue node node1 that is adjancent to node0, if node0 is the only white node adjacent to node1, then
     the forcing rule applies, and node0 will be forced blue in the next iteration, therefore return True. 
     If this doesn't apply for all blue nodes, then node0 won't be forced blue, therefore return False. """
    for node1 in set(graph.adj[node0]).intersection(b):
        if set(graph.adj[node1]).intersection(set(graph.nodes) - b) == {node0}:
            return True

    return False


def psd_rule(graph: nx.Graph, node0, b):
    """ Returns True if, given a graph and a set of blue vertices b, nodes0 will be forced blue in the next iteration
    using the PSD forcing rule.
    Returns False otherwise. """

    if graph.nodes[node0]["color"] == "blue":
        return True

    """ For every blue node node1 that is adjancent to node0, if node0 is not adjancent to any other white vertex 
    adjacent to node1, then node0 will be forced in the next iteration, therefore return True.
    If this doesn't apply for all blue nodes, then node0 won't be forced blue, therefore return False. """
    for node1 in set(graph.adj[node0]).intersection(b):
        works = True
        for node2 in set(graph.adj[node1]).intersection(set(graph.nodes) - b - {node0}):
            if node0 in set(graph.adj[node2]):
                works = False
                break

        if works:
            return True

    return False


def forcing_process(graph: nx.Graph, b, rule=forcing_rule, t=1):
    """ Applies the forcing rule given to the given graph with the given initial set of blue vertices b. Once the
    forcing process is done, returns the final set of blue vertices. """

    for node1 in b:
        graph.nodes[node1]["color"] = "blue"

    new_b = []
    for node1 in list(graph.nodes):
        if rule(graph, node1, set(b)):
            new_b.append(node1)

    if len(new_b) == len(list(graph.nodes)):
        return new_b, t

    if b == new_b:
        return new_b, t

    return forcing_process(graph, new_b, rule, t + 1)


def zero_forcing(graph: nx.Graph, rule=forcing_rule):
    """ Applies a Genetic Algorithm to the given graph in order to find its minimal zero-forcing set, and therefore
    zero-forcing number. Returns the zero-forcing number. (With a slight tweak the function can return the zero-forcing
    set aswell if need be)"""

    population_size = 8
    num_of_elite_chrom = 1
    tournament_selection_size = 4
    mutation_rate = 0.25
    target_gen = 500
    max_stagnant = 20 + len(list(graph.nodes)) * 2

    class Chromosome:
        def __init__(self):
            lst = list(graph.nodes)
            k = random.randrange(1, len(lst) + 1)
            self.genes = random.sample(lst, k)
            self.fitness = self.get_fitness()
            self.t = 0

        def get_genes(self):
            return self.genes

        def get_fitness(self):
            for node1 in graph.nodes:
                graph.nodes[node1]["color"] = "white"

            self.fitness = 0
            res = forcing_process(graph, self.genes, rule)
            self.t = res[1]
            forcing = len(res[0])
            if forcing != len(list(graph.nodes)):
                return -10 ** 20
            result = (10 ** 20) - (len(self.genes) + 2) ** 5 - self.t

            for node1 in graph.nodes:
                graph.nodes[node1]["color"] = "white"

            self.fitness = result
            return result

        def __str__(self):
            return self.genes.__str__()

    class Population:
        def __init__(self, size):
            self._chromosomes = []
            i = 0
            while i < size:
                self._chromosomes.append(Chromosome())
                i += 1

        def get_chromosomes(self): return self._chromosomes

    class GeneticAlgorithm:
        @staticmethod
        def evolve(pop):
            return GeneticAlgorithm._mutate_population(GeneticAlgorithm._crossover_population(pop))

        @staticmethod
        def _crossover_population(pop):
            crossover_pop = Population(0)
            for i in range(num_of_elite_chrom):
                crossover_pop.get_chromosomes().append(pop.get_chromosomes()[i])
            i = num_of_elite_chrom
            while i < population_size:
                chromosome1 = GeneticAlgorithm._select_tournament_population(pop).get_chromosomes()[0]
                chromosome2 = GeneticAlgorithm._select_tournament_population(pop).get_chromosomes()[0]
                crossover_pop.get_chromosomes().append(
                    GeneticAlgorithm._crossover_chromosomes(chromosome1, chromosome2))
                i += 1
            return crossover_pop

        @staticmethod
        def _mutate_population(pop):
            for i in range(num_of_elite_chrom, population_size):
                GeneticAlgorithm._mutate_chromosome(pop.get_chromosomes()[i])
            return pop

        @staticmethod
        def _crossover_chromosomes(chromosome1, chromosome2):
            crossover_chrom = Chromosome()
            k = random.randrange(min(len(chromosome1.genes), len(chromosome2.genes)), 1 + max(len(chromosome1.genes),
                                                                                              len(chromosome2.genes)))
            crossover_chrom.genes = list(set(random.sample(chromosome1.genes + chromosome2.genes, k)))
            crossover_chrom.fitness = crossover_chrom.get_fitness()

            return crossover_chrom

        @staticmethod
        def _mutate_chromosome(chromosome):
            if random.random() < mutation_rate:
                lst = list(graph.nodes) + chromosome.genes
                k = random.randrange(1, len(list(graph.nodes)) + 1)
                chromosome.genes = random.sample(lst, k)
                chromosome.genes = list(set(chromosome.genes))
                chromosome.fitness = chromosome.get_fitness()

        @staticmethod
        def _select_tournament_population(pop):
            tournament_pop = Population(0)
            i = 0
            while i < tournament_selection_size:
                tournament_pop.get_chromosomes().append(pop.get_chromosomes()[random.randrange(0, population_size)])
                i += 1
            tournament_pop.get_chromosomes().sort(key=lambda x: x.fitness, reverse=True)
            return tournament_pop

    def _print_population(pop, gen_number):
        print("\n--------------------------------------------------")
        print("Generation #", gen_number, "| Fittest chromosome fitness:", pop.get_chromosomes()[0].get_fitness())
        print("--------------------------------------------------")
        i = 0
        for x in pop.get_chromosomes():
            print("Chromosome  #", i, " :", x, "| Fitness: ", x.get_fitness())
            i += 1

    population = Population(population_size)
    population.get_chromosomes().sort(key=lambda x: x.get_fitness(), reverse=True)
    # _print_population(population, 0)
    generation_number = 1
    stagnant = 0
    while generation_number <= target_gen:
        old1, old2 = len(population.get_chromosomes()[0].genes), population.get_chromosomes()[0].t
        population = GeneticAlgorithm.evolve(population)
        population.get_chromosomes().sort(key=lambda x: x.get_fitness(), reverse=True)
        # _print_population(population, generation_number)
        if (len(population.get_chromosomes()[0].genes), population.get_chromosomes()[0].t) == (old1, old2):
            stagnant += 1
        else:
            stagnant = 0
        if stagnant > max_stagnant:
            population.get_chromosomes().sort(key=lambda x: x.fitness, reverse=True)
            return len(population.get_chromosomes()[0].genes), population.get_chromosomes()[0].t, \
                   population.get_chromosomes()[0].genes
        generation_number += 1

    population.get_chromosomes().sort(key=lambda x: x.fitness, reverse=True)
    return len(population.get_chromosomes()[0].genes), population.get_chromosomes()[0].t, population.get_chromosomes()[
        0].genes


def throttling_num(graph: nx.Graph, rule=forcing_rule):
    """ Applies a Genetic Algorithm to the given graph in order to find its minimal zero-forcing set, and therefore
    zero-forcing number. Returns the zero-forcing number. (With a slight tweak the function can return the zero-forcing
    set aswell if need be)"""

    population_size = 8
    num_of_elite_chrom = 1
    tournament_selection_size = 4
    mutation_rate = 0.25
    target_gen = 500
    max_stagnant = 20 + len(list(graph.nodes)) * 2

    class Chromosome:
        def __init__(self):
            lst = list(graph.nodes)
            k = random.randrange(1, len(lst) + 1)
            self.genes = random.sample(lst, k)
            self.fitness = self.get_fitness()
            self.t = 0

        def get_genes(self):
            return self.genes

        def get_fitness(self):
            for node1 in graph.nodes:
                graph.nodes[node1]["color"] = "white"

            self.fitness = 0
            res = forcing_process(graph, self.genes, rule)
            self.t = res[1]
            forcing = len(res[0])
            if forcing != len(list(graph.nodes)):
                return -10 ** 20
            result = (10 ** 20) - len(self.genes) - self.t

            for node1 in graph.nodes:
                graph.nodes[node1]["color"] = "white"

            self.fitness = result
            return result

        def __str__(self):
            return self.genes.__str__()

    class Population:
        def __init__(self, size):
            self._chromosomes = []
            i = 0
            while i < size:
                self._chromosomes.append(Chromosome())
                i += 1

        def get_chromosomes(self): return self._chromosomes

    class GeneticAlgorithm:
        @staticmethod
        def evolve(pop):
            return GeneticAlgorithm._mutate_population(GeneticAlgorithm._crossover_population(pop))

        @staticmethod
        def _crossover_population(pop):
            crossover_pop = Population(0)
            for i in range(num_of_elite_chrom):
                crossover_pop.get_chromosomes().append(pop.get_chromosomes()[i])
            i = num_of_elite_chrom
            while i < population_size:
                chromosome1 = GeneticAlgorithm._select_tournament_population(pop).get_chromosomes()[0]
                chromosome2 = GeneticAlgorithm._select_tournament_population(pop).get_chromosomes()[0]
                crossover_pop.get_chromosomes().append(
                    GeneticAlgorithm._crossover_chromosomes(chromosome1, chromosome2))
                i += 1
            return crossover_pop

        @staticmethod
        def _mutate_population(pop):
            for i in range(num_of_elite_chrom, population_size):
                GeneticAlgorithm._mutate_chromosome(pop.get_chromosomes()[i])
            return pop

        @staticmethod
        def _crossover_chromosomes(chromosome1, chromosome2):
            crossover_chrom = Chromosome()
            k = random.randrange(min(len(chromosome1.genes), len(chromosome2.genes)), 1 + max(len(chromosome1.genes),
                                                                                              len(chromosome2.genes)))
            crossover_chrom.genes = list(set(random.sample(chromosome1.genes + chromosome2.genes, k)))
            crossover_chrom.fitness = crossover_chrom.get_fitness()

            return crossover_chrom

        @staticmethod
        def _mutate_chromosome(chromosome):
            if random.random() < mutation_rate:
                lst = list(graph.nodes) + chromosome.genes
                k = random.randrange(1, len(list(graph.nodes)) + 1)
                chromosome.genes = random.sample(lst, k)
                chromosome.genes = list(set(chromosome.genes))
                chromosome.fitness = chromosome.get_fitness()

        @staticmethod
        def _select_tournament_population(pop):
            tournament_pop = Population(0)
            i = 0
            while i < tournament_selection_size:
                tournament_pop.get_chromosomes().append(pop.get_chromosomes()[random.randrange(0, population_size)])
                i += 1
            tournament_pop.get_chromosomes().sort(key=lambda x: x.fitness, reverse=True)
            return tournament_pop

    def _print_population(pop, gen_number):
        print("\n--------------------------------------------------")
        print("Generation #", gen_number, "| Fittest chromosome fitness:", pop.get_chromosomes()[0].get_fitness())
        print("--------------------------------------------------")
        i = 0
        for x in pop.get_chromosomes():
            print("Chromosome  #", i, " :", x, "| Fitness: ", x.get_fitness())
            i += 1

    population = Population(population_size)
    population.get_chromosomes().sort(key=lambda x: x.get_fitness(), reverse=True)
    # _print_population(population, 0)
    generation_number = 1
    stagnant = 0
    while generation_number <= target_gen:
        old1, old2 = len(population.get_chromosomes()[0].genes), population.get_chromosomes()[0].t
        population = GeneticAlgorithm.evolve(population)
        population.get_chromosomes().sort(key=lambda x: x.get_fitness(), reverse=True)
        # _print_population(population, generation_number)
        if (len(population.get_chromosomes()[0].genes), population.get_chromosomes()[0].t) == (old1, old2):
            stagnant += 1
        else:
            stagnant = 0
        if stagnant > max_stagnant:
            population.get_chromosomes().sort(key=lambda x: x.fitness, reverse=True)
            return len(population.get_chromosomes()[0].genes) + population.get_chromosomes()[0].t
        generation_number += 1

    population.get_chromosomes().sort(key=lambda x: x.fitness, reverse=True)
    return len(population.get_chromosomes()[0].genes) + population.get_chromosomes()[0].t


G = nx.Graph()
G.add_edges_from([(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7)])
for node in G.nodes:
    G.nodes[node]["color"] = "white"

z = zero_forcing(G, psd_rule)
t = throttling_num(G, psd_rule)
print("1. Zero Frocing number: " + str(z[0]) + "  |  Propagation number: " + str(z[1]) + "  |  Zero Forcing Set: " +
      str(z[2]))
print("2. Minimal throttling number: " + str(t))


"""A = GraphMatrix(G, [1, 1, 1, 1, 1, 1, 1])
A.show(True)
print("\n\nMatrix Rank: " + str(A.rank()))
"""
