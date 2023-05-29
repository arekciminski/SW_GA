import pygad

from Get_data_from_DWDS import Epa
import numpy as np

class GA_function:

    def __init__(self, ga_parameters, sw_parameters):
        self.sw_parameters = sw_parameters
        self.ga_parameters = ga_parameters
        self.ep = Epa(sw_parameters)

    def fitnes_function(self, ga_instance, species, solution_idx):
        self.ep.get_data(species)

        fitnes_fun = self.ga_parameters['penalty_function_weights'][0] *self.error_penalty_function(self.ep) + \
            self.ga_parameters['penalty_function_weights'][0] * self.head_penalty_function(self.ep) + \
            self.ga_parameters['penalty_function_weights'][0] * self.flow_penalty_function(self.ep) + \
            self.ga_parameters['penalty_function_weights'][0] * self.tank_penalty_function(self.ep) + \
            self.ga_parameters['penalty_function_weights'][0] * self.tank_final_state_penalty_function(self.ep)
        return -fitnes_fun

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

        for j in range(epa.sw_parameters['time_duration_h']):
            for i in range(len(epa.sw_parameters['mes_nodes_names']) - len(epa.sw_parameters['tanks_names'])):
                cross_values = -epa.data['head_output_' + epa.sw_parameters['mes_nodes_names'][i]][0][j] + \
                               epa.sw_parameters['min_head_lev'][i]
                if cross_values < 0:
                    head_penalty += abs(cross_values)

                cross_values = epa.data['head_output_' + epa.sw_parameters['mes_nodes_names'][i]][0][j] - \
                               epa.sw_parameters['max_head_lev'][i]
                if cross_values < 0:
                    head_penalty += abs(cross_values)

        return head_penalty

    def tank_penalty_function(self, epa):
        tank_penalty = 0

        for j in range(epa.sw_parameters['time_duration_h']):
            for i in range(len(epa.sw_parameters['mes_nodes_names']) - len(epa.sw_parameters['tanks_names']),
                           len(epa.sw_parameters['mes_nodes_names'])):
                cross_values = -epa.data['head_output_' + epa.sw_parameters['mes_nodes_names'][i]][0][j] + \
                               epa.sw_parameters['min_head_lev'][i]
                if cross_values < 0:
                    tank_penalty += abs(cross_values)

                cross_values = epa.data['head_output_' + epa.sw_parameters['mes_nodes_names'][i]][0][j] - \
                               epa.sw_parameters['max_head_lev'][i]
                if cross_values < 0:
                    tank_penalty += abs(cross_values)

        return tank_penalty

    def flow_penalty_function(self, epa):

        flow_penalty = 0
        for j in range(epa.sw_parameters['time_duration_h']):
            for i in range(len(epa.sw_parameters['mes_links_names'])):
                cross_values = -epa.data['flow_output_' + epa.sw_parameters['mes_links_names'][i]][0][j] + \
                               epa.sw_parameters['min_flow_lev'][i]
                if cross_values < 0:
                    flow_penalty += abs(cross_values)

            for i in range(len(epa.sw_parameters['mes_links_names'])):
                cross_values = epa.data['flow_output_' + epa.sw_parameters['mes_links_names'][i]][0][j] - \
                               epa.sw_parameters['max_flow_lev'][i]
                if cross_values < 0:
                    flow_penalty += abs(cross_values)
        return flow_penalty

    def error_penalty_function(self,epa):
        return sum(epa.data['error_output'][0])

    def GA_mutation(self,pop):

        low_int_val = self.ga_parameters['pop_int_range_low']
        high_int_val = self.ga_parameters['pop_int_range_high']
        low_float_val = self.ga_parameters['pop_float_range_low']
        high_float_val = self.ga_parameters['pop_float_range_high']

        pop['mate_mut'] = pop['pop']
        if self.ga_parameters['mutation_type'] == 'random':
            mate =[]

            for i in range(self.ga_parameters['num_specimen']):

                if self.ga_parameters['num_int_genes'] > 0:
                    if np.random.rand() < self.ga_parameters['mutation_probability']/100:
                        num_genes = int(self.ga_parameters['num_int_genes']*\
                                        self.ga_parameters['mutation_percent_genes']/100)
                        genes_position = np.random.randint(self.ga_parameters['num_int_genes'],size=(num_genes,))
                        for ge_pos in genes_position:
                            pop['mate_mut'][i][0][ge_pos] = np.random.randint(low_int_val, high_int_val,size=(1,))

                if self.ga_parameters['num_float_genes'] > 0:
                    if np.random.rand() < self.ga_parameters['mutation_probability']/100:
                       num_genes = int(
                            self.ga_parameters['num_int_genes'] * self.ga_parameters['mutation_percent_genes'] / 100)
                       genes_position = np.random.randint(self.ga_parameters['num_int_genes'], size=(num_genes,))
                       for ge_pos in genes_position:
                           pop['mate_mut'][i][1][ge_pos] = (high_float_val - low_float_val) * np.random.random() + \
                                                           low_float_val
        return 0

    def GA_crossover(self, pop):

        return 0

    def GA_selection(self,pop):

        return 0

    def GA_initialization(self):

        low_int_val = self.ga_parameters['pop_int_range_low']
        high_int_val = self.ga_parameters['pop_int_range_high']
        low_float_val = self.ga_parameters['pop_float_range_low']
        high_float_val = self.ga_parameters['pop_float_range_high']

        pop = {'pop' : [],
               'fitnes_function' : []}

        for i in range(self.ga_parameters['num_specimen']):
            specimen = []

            if self.ga_parameters['num_int_genes']>0:
                specimen.append((high_int_val - low_int_val) * np.random.randint(low_int_val , high_int_val ,size = (self.ga_parameters['num_int_genes'],)) + low_int_val)

            if self.ga_parameters['num_float_genes']>0:
                specimen.append((high_float_val - low_float_val) * np.random.random_sample(\
                    (self.ga_parameters['num_float_genes'],)) + low_float_val)

            pop['pop'].append(specimen)
            pop['fitnes_function'].append(0)

        return pop


if __name__ == '__main__':
        print('\nTo jest biblioteka pomocna przy algorytmie genetycznym')


        file_name = "chojnice_kwiecien_obiekt.inp"

        sw_parameters = {
            'file_name': file_name,
            'mes_links_names': ['F1', 'K1', 'P1'],
            'mes_nodes_names': ['001', '002', '055', '083', '158', '180'],
            'pumps_names': ['F1', 'K1', 'P1'],
            'pumps_patterns_names': ['f1', 'k1', 'p1'],
            'demand_patterns_names': ['par2', 'kar2'],
            'time_duration_h': 24,  # [h]
            'hydraulic_step_s': 3600,  # [h]
            'tanks_names': ['180'],
            'initial_tanks_lev': [3.58],
            'hydraulic_monitor_values': ['flow', 'energy', 'head'],
            'min_head_lev': [0, 0, 0, 199, 195, 167.2],
            'max_head_lev': [200, 200, 200, 210, 210, 170.8],
            'min_flow_lev': [0, 0, 0],
            'max_flow_lev': [330, 400, 140],

        }

        sw_parameters['num_mes_links'] = len(sw_parameters['mes_links_names'])
        sw_parameters['num_mes_nodes'] = len(sw_parameters['mes_nodes_names'])
        sw_parameters['num_pumps'] = len(sw_parameters['pumps_names'])
        sw_parameters['num_demands'] = len(sw_parameters['demand_patterns_names'])
        sw_parameters['num_tanks'] = len(sw_parameters['tanks_names'])
        sw_parameters['time_duration_s'] = sw_parameters['time_duration_h'] * 3600  ##[s]

        ga_parameters = {
            'penalty_function_weights': [1e5, 10, 20, 30, 40],
            'num_generations': 100,
            'num_specimen': 10,
            'gene_type': float,
            'num_float_genes': sw_parameters['num_pumps'] * sw_parameters['time_duration_h'],
            'pop_int_range_low': 0,
            'pop_int_range_high': 1,
            'num_int_genes': 10,
            'parent_selection_type': 'sss',  # rws, sus, random, tournament
            'crossover_probability': 0.9,
            'crossover_type': 'single_point',  # two_points , uniform , scattered
            'mutation_type': 'random',  # swap, inversion , scramble ,
            'mutation_percent_genes': 50,
            'mutation_probability': 70,
            'pop_float_range_low': 0,
            'pop_float_range_high': 1,
        }

        ga = GA_function(ga_parameters , sw_parameters)

        popul = ga.GA_initialization()

        popul  = ga.GA_mutation(popul)

