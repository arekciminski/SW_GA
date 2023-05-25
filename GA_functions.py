import pygad
from Get_data_from_DWDS import Epa

class GA_function:

    def __init__(self, GA_parameters, SW_parameters):
        self.SW_parameters = SW_parameters
        self.GA_parameters = GA_parameters
        self.ep = Epa(SW_parameters)

    def fitnes_function(self, species):
        self.ep.get_data(species)
        print(self.error_penalty_function(self.ep))
        print(self.head_penalty_function(self.ep))
        print(self.flow_penalty_function(self.ep))
        print(self.tank_penalty_function(self.ep))
        print(self.tank_final_state_penalty_function(self.ep))
        print(self.ep.data.keys())

        return 0

    def GA_loop(self):
        ga_instance = pygad.GA(num_generations=self.GA_parameters['num_generations'],
                           num_parents_mating=self.GA_parameters['num_parents_mating'],
                           fitness_func = self.fitnes_function(),
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
                           random_mutation_min_val = self.GA_parameters['random_mutation_min_val'],
                           random_mutation_max_val = self.GA_parameters['random_mutation_max_val'],
                           on_mutation = self.GA_parameters['on_mutation'],
                           save_best_solutions = self.GA_parameters['save_best_solutions'],
                           parallel_processing = self.GA_parameters['parallel_processing'],
        )

    def tank_final_state_penalty_function(self, epa):
        for i in range(len(epa.SW_parameters['mes_nodes_names']) - len(epa.SW_parameters['tanks_names']),
                       len(epa.SW_parameters['mes_nodes_names'])):
            tank_inital_state = epa.data['head_output_' + epa.SW_parameters['mes_nodes_names'][i]][0][0]
            tank_final_state = epa.data['head_output_' + epa.SW_parameters['mes_nodes_names'][i]][0][-1]

        return abs(tank_inital_state - tank_final_state)

    def tank_final_state_penalty_function(self, epa):

        for i in range(len(epa.SW_parameters['mes_nodes_names']) - len(epa.SW_parameters['tanks_names']),
                       len(epa.SW_parameters['mes_nodes_names'])):
            tank_inital_state = epa.data['head_output_' + epa.SW_parameters['mes_nodes_names'][i]][0][0]
            tank_final_state = epa.data['head_output_' + epa.SW_parameters['mes_nodes_names'][i]][0][-1]

        return abs(tank_inital_state - tank_final_state)

    def head_penalty_function(self, epa):

        head_penalty = 0

        for j in range(epa.SW_parameters['time_duration_h']):
            for i in range(len(epa.SW_parameters['mes_nodes_names']) - len(epa.SW_parameters['tanks_names'])):
                cross_values = -epa.data['head_output_' + epa.SW_parameters['mes_nodes_names'][i]][0][j] + \
                               epa.SW_parameters['min_head_lev'][i]
                if cross_values < 0:
                    head_penalty += abs(cross_values)

                cross_values = epa.data['head_output_' + epa.SW_parameters['mes_nodes_names'][i]][0][j] - \
                               epa.SW_parameters['max_head_lev'][i]
                if cross_values < 0:
                    head_penalty += abs(cross_values)

        return head_penalty

    def tank_penalty_function(self, epa):
        tank_penalty = 0

        for j in range(epa.SW_parameters['time_duration_h']):
            for i in range(len(epa.SW_parameters['mes_nodes_names']) - len(epa.SW_parameters['tanks_names']),
                           len(epa.SW_parameters['mes_nodes_names'])):
                cross_values = -epa.data['head_output_' + epa.SW_parameters['mes_nodes_names'][i]][0][j] + \
                               epa.SW_parameters['min_head_lev'][i]
                if cross_values < 0:
                    tank_penalty += abs(cross_values)

                cross_values = epa.data['head_output_' + epa.SW_parameters['mes_nodes_names'][i]][0][j] - \
                               epa.SW_parameters['max_head_lev'][i]
                if cross_values < 0:
                    tank_penalty += abs(cross_values)

        return tank_penalty

    def flow_penalty_function(self, epa):

        flow_penalty = 0
        for j in range(epa.SW_parameters['time_duration_h']):
            for i in range(len(epa.SW_parameters['mes_links_names'])):
                cross_values = -epa.data['flow_output_' + epa.SW_parameters['mes_links_names'][i]][0][j] + \
                               epa.SW_parameters['min_flow_lev'][i]
                if cross_values < 0:
                    flow_penalty += abs(cross_values)

            for i in range(len(epa.SW_parameters['mes_links_names'])):
                cross_values = epa.data['flow_output_' + epa.SW_parameters['mes_links_names'][i]][0][j] - \
                               epa.SW_parameters['max_flow_lev'][i]
                if cross_values < 0:
                    flow_penalty += abs(cross_values)
        return flow_penalty

    def error_penalty_function(self,epa):
        return sum(epa.data['error_output'][0])


if __name__ == '__main__':
        print('To jest biblioteka pomocna przy algorytmie genetycznym')