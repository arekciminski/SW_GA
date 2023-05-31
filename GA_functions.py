from Get_data_from_DWDS import Epa
import numpy as np
import copy
from statistics import mean
from tqdm import tqdm


class GA_function:

    def __init__(self, ga_parameters, sw_parameters):
        self.sw_parameters = sw_parameters
        self.ga_parameters = ga_parameters
        self.ep = Epa(sw_parameters)
        self.pop = []

    def tank_final_state_penalty_function(self, epa):
        for i in range(len(epa.sw_parameters['mes_nodes_names']) - len(epa.sw_parameters['tanks_names']),
                       len(epa.sw_parameters['mes_nodes_names'])):
            tank_inital_state = epa.data['head_output_' + epa.sw_parameters['mes_nodes_names'][i]][0][0]
            tank_final_state = epa.data['head_output_' + epa.sw_parameters['mes_nodes_names'][i]][0][-1]

        return abs(tank_inital_state - tank_final_state)

    def tank_final_state_penalty_function(self, epa):

        for i in range(len(epa.sw_parameters['mes_nodes_names']) - len(epa.sw_parameters['tanks_names']),
                       len(epa.sw_parameters['mes_nodes_names'])):
            tank_inital_state = epa.data['head_output_' + epa.sw_parameters['mes_nodes_names'][i]][0][0]
            tank_final_state = epa.data['head_output_' + epa.sw_parameters['mes_nodes_names'][i]][0][-1]

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

        error_values = epa.data['error_output'][0]
        for i in range(len(error_values)):
            if error_values[i] > 0:
                error_values[i] = 1

        return sum(epa.data['error_output'][0])

    def ga_mutation(self):

        low_int_val = self.ga_parameters['pop_int_range_low']
        high_int_val = self.ga_parameters['pop_int_range_high']
        low_float_val = self.ga_parameters['pop_float_range_low']
        high_float_val = self.ga_parameters['pop_float_range_high']
        self.pop['mate_mut'] = []
        specimen = copy.deepcopy(self.pop['pop'])
        if self.ga_parameters['mutation_type'] == 'random':
            for i in range(self.ga_parameters['num_specimen']):
                if np.random.rand() < self.ga_parameters['mutation_percent_probability'] / 100:
                    species = specimen[i]
                    if self.ga_parameters['num_int_genes'] > 0:

                            num_genes = int(self.ga_parameters['num_int_genes']*\
                                            self.ga_parameters['mutation_percent_genes']/100)

                            genes_position = np.random.randint(self.ga_parameters['num_int_genes'],size=(num_genes,))

                            for ge_pos in genes_position:
                                species[1][ge_pos] = np.random.randint(low_int_val, high_int_val+1,size=(1,))
                    if self.ga_parameters['num_float_genes'] > 0:
                        if np.random.rand() < self.ga_parameters['mutation_percent_probability']/100:

                           num_genes = int(
                                self.ga_parameters['num_float_genes'] * self.ga_parameters['mutation_percent_genes'] / 100)

                           genes_position = np.random.randint(self.ga_parameters['num_float_genes'], size=(num_genes,))

                           for ge_pos in genes_position:
                               species[2][ge_pos] = (high_float_val - low_float_val) * np.random.random() + low_float_val
                    self.pop['mate_mut'].append(species)
        return 1

    def ga_crossover(self):

        #low_int_val = self.ga_parameters['pop_int_range_low']
        #high_int_val = self.ga_parameters['pop_int_range_high']
        #low_float_val = self.ga_parameters['pop_float_range_low']
        #high_float_val = self.ga_parameters['pop_float_range_high']

        self.pop['mate_cros'] = []
        if self.ga_parameters['crossover_type'] == 'single_point':
            for i in range(self.ga_parameters['num_specimen']):
                specimen_to_mutation = np.random.randint(1,self.ga_parameters['num_specimen'],size=(2,))

                specimen_0 = list(copy.deepcopy(self.pop['pop'][specimen_to_mutation[0]]))
                specimen_1 = list(copy.deepcopy(self.pop['pop'][specimen_to_mutation[1]]))

                if self.ga_parameters['num_int_genes'] > 0:

                    cros_position = np.random.randint(1,self.ga_parameters['num_int_genes'], size=(1,))

                    first_species_first_half = specimen_0[1][:cros_position[0]]
                    first_species_second_half = specimen_0[1][cros_position[0]:]
                    second_species_first_half = specimen_1[1][:cros_position[0]]
                    second_species_second_half = specimen_1[1][cros_position[0]:]

                    specimen_0[1] = np.concatenate((second_species_first_half, first_species_second_half))
                    specimen_1[1] = np.concatenate((first_species_first_half, second_species_second_half))

                if self.ga_parameters['num_float_genes'] > 0:

                    cros_position = np.random.randint(1 , self.ga_parameters['num_float_genes'], size=(1,))

                    first_species_first_half = specimen_0[2][:cros_position[0]]
                    first_species_second_half = specimen_0[2][cros_position[0]:]
                    second_species_first_half = specimen_1[2][:cros_position[0]]
                    second_species_second_half = specimen_1[2][cros_position[0]:]

                    specimen_0[2] = np.concatenate((second_species_first_half , first_species_second_half))
                    specimen_1[2] = np.concatenate((first_species_first_half, second_species_second_half))

                self.pop['mate_cros'].append(specimen_0)
                self.pop['mate_cros'].append(specimen_1)

        return 1

    def ga_selection(self):

        if self.ga_parameters['parent_selection_type'] == 'first_n':
            keys_list = ['pop', 'mate_mut', 'mate_cros']
            tuple_list = []
            for key in keys_list:
                for i in range(len(self.pop[key])):
                    tuple_list.append(tuple((self.pop[key][i][0],self.pop[key][i][1],self.pop[key][i][2])))

            sorted_fitnes_fun = sorted(tuple_list, key=lambda x: x[0])[:self.ga_parameters['num_specimen']]

            self.pop['pop'] = []

            for i in range(self.ga_parameters['num_specimen']):
                temp = []
                for j in range(len(sorted_fitnes_fun[0])):
                    temp.append(sorted_fitnes_fun[i][j])
                self.pop['pop'].append(temp)
            self.ga_statistics()

        return 0

    def ga_statistics(self):

        self.pop['num_iteration'] += 1

        self.pop['last_best_solution'] = self.pop['pop'][0]

        self.pop['best_solutions'].append(self.pop['last_best_solution'][0])

        population_fun_fitnes = []

        for i in range(ga_parameters['num_specimen']):
            population_fun_fitnes.append(self.pop['pop'][i][0])

        self.pop['mean_pop_value'].append(mean(population_fun_fitnes))

    def ga_initialization(self):

        low_int_val = self.ga_parameters['pop_int_range_low']
        high_int_val = self.ga_parameters['pop_int_range_high']
        low_float_val = self.ga_parameters['pop_float_range_low']
        high_float_val = self.ga_parameters['pop_float_range_high']

        self.pop = {'pop': []}
        self.pop['num_iteration'] = 0
        self.pop['best_solutions'] = []
        self.pop['mean_pop_value'] = []

        for i in range(self.ga_parameters['num_specimen']):
            specimen = [[], [], []]

            if self.ga_parameters['num_int_genes'] > 0:
                specimen[1] = (high_int_val - low_int_val) * np.random.randint(low_int_val, high_int_val + 1, size= \
                    (self.ga_parameters['num_int_genes'],)) + low_int_val
            if self.ga_parameters['num_float_genes'] > 0:
                specimen[2] = (high_float_val - low_float_val) * np.random.random_sample(
                    (self.ga_parameters['num_float_genes'],)) + low_float_val
            self.pop['pop'].append(specimen)

        return self.pop

    def ga_fitnes_function(self):

        if self.pop['num_iteration'] == 0:
            keys_list = ['pop']
        else:
            keys_list = ['pop', 'mate_mut', 'mate_cros']
        for key in keys_list:
            for i in tqdm(range(len(self.pop[key]))):
                species = self.pop[key][i][2]
                self.ep.get_data(species)

                fitnes_fun = self.ga_parameters['penalty_function_weights'][0] * self.error_penalty_function(self.ep) +\
                    self.ga_parameters['penalty_function_weights'][1] * self.head_penalty_function(self.ep) +\
                    self.ga_parameters['penalty_function_weights'][2] * self.flow_penalty_function(self.ep) +\
                    self.ga_parameters['penalty_function_weights'][3] * self.tank_penalty_function(self.ep) +\
                    self.ga_parameters['penalty_function_weights'][4] * self.tank_final_state_penalty_function(self.ep)

                self.pop[key][i][0] = fitnes_fun
        return self.pop

    def ga_stop_criterion(self):

        if self.pop['num_iteration'] > self.ga_parameters['stop_criterion_min_generations']:

            diff_fitness_function = []

            last_solutions = self.pop['best_solutions'][-self.ga_parameters['stop_criterion_min_generations']:]

            for i in range(self.ga_parameters['stop_criterion_min_generations']-1, 0,-1):

                if abs(last_solutions[i-1] - last_solutions[i]) < self.ga_parameters['stop_criterion_min_change']:
                    diff_fitness_function.append(0)
                else:
                    diff_fitness_function.append(1)

            return sum(diff_fitness_function)


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
            'penalty_function_weights': [1e9, 1, 1, 1, 1],
            'num_generations': 200,
            'num_specimen': 30,
            'num_float_genes': sw_parameters['num_pumps'] * sw_parameters['time_duration_h'],
            'pop_int_range_low': 0,
            'pop_int_range_high': 1,
            'num_int_genes': 10,
            'parent_selection_type': 'first_n',  # rws, sus, random, tournament
            'crossover_percent_probability': 40,
            'crossover_type': 'single_point',  # two_points , uniform , scattered
            'mutation_type': 'random',  # swap, inversion , scramble ,
            'mutation_percent_genes': 40,
            'mutation_percent_probability': 80,
            'pop_float_range_low': 0.4,
            'pop_float_range_high': 1,
            'stop_criterion_min_generations': 15,
            'stop_criterion_min_change': 1e-3,
        }

        ga = GA_function(ga_parameters , sw_parameters)

        ga.ga_initialization()

        ga.ga_fitnes_function()

        for i in range(ga.ga_parameters['num_generations']):

            ga.ga_mutation()

            ga.ga_crossover()

            ga.ga_fitnes_function()

            ga.ga_selection()

            print(f"Iteration number: {ga.pop['num_iteration']}, Mean population value: {ga.pop['mean_pop_value'][-1]}, Best solution: {ga.pop['best_solutions'][-1]}")

            if ga.ga_stop_criterion() == 0:
                break


        dump