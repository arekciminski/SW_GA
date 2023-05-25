import pygad

class GA_functions:

    def __init__(self, species, GA_parameters, SW_parameters):
        self.SW_parameters = SW_parameters
        self.GA_parameters = GA_parameters
        self.species = species

    def fitnes_function(self):

        return 0

    def GA_loop(self):
        ga_instance = pygad.GA(num_generations=self.GA_parameters['num_generations'],
                           num_parents_mating=self.GA_parameters['num_parents_mating'],
                           fitness_func = fitnes_function(),
                           sol_per_pop=self.GA_parameters['sol_per_pop'],
                           num_genes=self.GA_parameters['num_genes'],
                           init_range_low=self.GA_parameters['init_range_low'],
                           init_range_high=self.GA_parameters['init_range_high'],
                           parent_selection_type=self.GA_parameters['parent_selection_type'],
                           keep_parents=self.GA_parameters['keep_parents'],
                           crossover_type=self.GA_parameters['crossover_type'],
                           mutation_type=self.GA_parameters['mutation_type'],
                           mutation_percent_genes=self.GA_parameters['mutation_percent_genes)'],
                           keep_elitism  =  self.GA_parameters['keep_elitism'],
                           crossover_probability = self.GA_parameters['crossover_probability'],
                           mutation_percent_genes = self.GA_parameters['mutation_percent_genes'],
                           random_mutation_min_val = self.GA_parameters['random_mutation_min_val'],
                           random_mutation_max_val = self.GA_parameters['random_mutation_max_val'],
                           on_mutation = self.GA_parameters['on_mutation'],
                           save_best_solutions = self.GA_parameters['save_best_solutions'],
                           parallel_processing = self.GA_parameters['parallel_processing'],
        )