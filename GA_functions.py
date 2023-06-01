import copy
import random
from statistics import mean

import numpy as np
from tqdm import tqdm

from Get_data_from_DWDS import Epa


class GA_function:

    def __init__(self, ga_parameters, sw_parameters):
        self.sw_parameters = sw_parameters
        self.ga_parameters = ga_parameters
        self.ep = Epa(sw_parameters)
        self.pop = []

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

        '''        for i in range(len(error_values)):
            if error_values[i] > 0:
                error_values[i] = 1'''

        return sum(error_values)

    def ga_mutation(self):

        rand = random.random()

        if rand < 0.5:
            if self.ga_parameters['mutation_type'] == 'random':
                self.ga_mutation_random()

        else:
            self.ga_mutiation_SGO_error()
            self.ga_mutiation_SGO_head()

    def ga_mutation_random(self):

        low_int_val = self.ga_parameters['pop_int_range_low']
        high_int_val = self.ga_parameters['pop_int_range_high']
        low_float_val = self.ga_parameters['pop_float_range_low']
        high_float_val = self.ga_parameters['pop_float_range_high']

        self.pop['mate_mut'] = []
        specimen = copy.deepcopy(self.pop['pop'])

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

    def ga_mutiation_SGO_error(self):
        low_int_val = self.ga_parameters['pop_int_range_low']
        high_int_val = self.ga_parameters['pop_int_range_high']
        low_float_val = self.ga_parameters['pop_float_range_low']
        high_float_val = self.ga_parameters['pop_float_range_high']

        self.pop['mate_mut'] = []
        specimen = copy.deepcopy(self.pop['pop'])


        time_dur = self.sw_parameters['time_duration_h']
        for i in range(self.ga_parameters['num_specimen']):
            if np.random.rand() < self.ga_parameters['mutation_percent_probability'] / 100:
                species = specimen[i]
                if self.ga_parameters['num_float_genes'] > 0:
                    for j in range(len(species[3])):
                        error_value = species[3][j]
                        flow_pump_karolewo = species[6][0][j]
                        flow_pump_funka = species[6][1][j]
                        flow_pump_plac = species[6][2][j]
                        if error_value == 4 and flow_pump_karolewo == 0 and flow_pump_plac > 0:
                            species[2][j + time_dur] += np.random.rand() * self.ga_parameters[
                                'specialize_operators_error_delta_value']
                            species[2][j + 2 * time_dur] -= np.random.rand() * self.ga_parameters[
                                'specialize_operators_error_delta_value']
                        if error_value == 4 and flow_pump_funka == 0:
                            species[2][j] += np.random.rand() * self.ga_parameters[
                                'specialize_operators_error_delta_value']
                        if error_value == 4 and flow_pump_plac == 0 and flow_pump_karolewo > 0:
                            species[2][j + 2 * time_dur] += np.random.rand() * self.ga_parameters[
                                'specialize_operators_error_delta_value']
                            species[2][j + time_dur] -= np.random.rand() * self.ga_parameters[
                                'specialize_operators_error_delta_value']
                        if error_value == 4 and flow_pump_plac == 0 and flow_pump_karolewo == 0:
                            species[2][j + 2 * time_dur] += np.random.rand() * self.ga_parameters[
                                'specialize_operators_error_delta_value']
                            species[2][j + time_dur] += np.random.rand() * self.ga_parameters[
                                'specialize_operators_error_delta_value']

                        if species[2][j] > high_float_val:
                            species[2][j] = high_float_val
                        if species[2][j] < low_float_val:
                            species[2][j] = low_float_val

                self.pop['mate_mut'].append(species)

    def ga_mutiation_SGO_head(self):
        low_int_val = self.ga_parameters['pop_int_range_low']
        high_int_val = self.ga_parameters['pop_int_range_high']
        low_float_val = self.ga_parameters['pop_float_range_low']
        high_float_val = self.ga_parameters['pop_float_range_high']

        self.pop['mate_mut'] = []
        specimen = copy.deepcopy(self.pop['pop'])


        time_dur = self.sw_parameters['time_duration_h']
        for i in range(self.ga_parameters['num_specimen']):
            if np.random.rand() < self.ga_parameters['mutation_percent_probability'] / 100:
                species = specimen[i]
                if self.ga_parameters['num_float_genes'] > 0:
                    for j in range(len(species[3])):
                        for k in range(len(self.sw_parameters['mes_nodes_names']) - \
                                       len(self.sw_parameters['tanks_names'])):
                            head_WMO = species[4][k][j]
                            head_min_lev = self.sw_parameters['min_head_lev'][k]
                            head_max_level = self.sw_parameters['max_head_lev'][k]
                            delta_min_head = head_min_lev - head_WMO
                            delta_max_head = head_WMO - head_max_level
                            if k == 0:
                                if  delta_min_head < 0:
                                    species[2][j] += np.random.rand() * self.ga_parameters[
                                    'specialize_operators_head_delta_value']
                                if delta_max_head < 0:
                                    species[2][j] -= np.random.rand() * self.ga_parameters[
                                        'specialize_operators_head_delta_value']
                            if k == 1:
                                if  delta_min_head < 0:
                                    species[2][j + time_dur] += np.random.rand() * self.ga_parameters[
                                    'specialize_operators_head_delta_value']
                                if delta_max_head < 0:
                                    species[2][j + time_dur] -= np.random.rand() * self.ga_parameters[
                                        'specialize_operators_head_delta_value']

                            if k == 2:
                                if  delta_min_head < 0:
                                    species[2][j + 2 * time_dur] += np.random.rand() * self.ga_parameters[
                                    'specialize_operators_head_delta_value']
                                if delta_max_head < 0:
                                    species[2][j + 2 * time_dur] -= np.random.rand() * self.ga_parameters[
                                        'specialize_operators_head_delta_value']

                            if k == 3:
                                if  delta_min_head < 0:
                                    species[2][j + time_dur] += np.random.rand() * self.ga_parameters[
                                    'specialize_operators_head_delta_value']/2
                                    species[2][j + 2 * time_dur] += np.random.rand() * self.ga_parameters[
                                    'specialize_operators_head_delta_value']
                                if delta_max_head < 0:
                                    species[2][j + time_dur] -= np.random.rand() * self.ga_parameters[
                                        'specialize_operators_head_delta_value']/2
                                    species[2][j + 2 * time_dur] -= np.random.rand() * self.ga_parameters[
                                        'specialize_operators_head_delta_value']

                            if k == 4:
                                if  delta_min_head < 0:
                                    species[2][j + time_dur] += np.random.rand() * self.ga_parameters[
                                    'specialize_operators_head_delta_value']
                                    species[2][j + 2 * time_dur] += np.random.rand() * self.ga_parameters[
                                    'specialize_operators_head_delta_value']/2
                                if delta_max_head < 0:
                                    species[2][j + time_dur] -= np.random.rand() * self.ga_parameters[
                                        'specialize_operators_head_delta_value']
                                    species[2][j + 2 * time_dur] -= np.random.rand() * self.ga_parameters[
                                        'specialize_operators_head_delta_value']/2


                        if species[2][j] > high_float_val:
                            species[2][j] = high_float_val
                        if species[2][j] < low_float_val:
                            species[2][j] = low_float_val

                self.pop['mate_mut'].append(species)


    def ga_crossover(self):

        #low_int_val = self.ga_parameters['pop_int_range_low']
        #high_int_val = self.ga_parameters['pop_int_range_high']
        #low_float_val = self.ga_parameters['pop_float_range_low']
        #high_float_val = self.ga_parameters['pop_float_range_high']

        self.pop['mate_cros'] = []
        if self.ga_parameters['crossover_type'] == 'single_point':
            for i in range(self.ga_parameters['num_specimen']):
                if np.random.rand() < self.ga_parameters['crossover_percent_probability'] / 100:
                    specimen_to_mutation = np.random.randint(1, self.ga_parameters['num_specimen'], size=(2,))

                    specimen_0 = list(copy.deepcopy(self.pop['pop'][specimen_to_mutation[0]]))
                    specimen_1 = list(copy.deepcopy(self.pop['pop'][specimen_to_mutation[1]]))

                    if self.ga_parameters['num_int_genes'] > 0:
                        cros_position = np.random.randint(1, self.ga_parameters['num_int_genes'], size=(1,))

                        first_species_first_half = specimen_0[1][:cros_position[0]]
                        first_species_second_half = specimen_0[1][cros_position[0]:]
                        second_species_first_half = specimen_1[1][:cros_position[0]]
                        second_species_second_half = specimen_1[1][cros_position[0]:]

                        specimen_0[1] = np.concatenate((second_species_first_half, first_species_second_half))
                        specimen_1[1] = np.concatenate((first_species_first_half, second_species_second_half))

                    if self.ga_parameters['num_float_genes'] > 0:
                        cros_position = np.random.randint(1, self.ga_parameters['num_float_genes'], size=(1,))

                        first_species_first_half = specimen_0[2][:cros_position[0]]
                        first_species_second_half = specimen_0[2][cros_position[0]:]
                        second_species_first_half = specimen_1[2][:cros_position[0]]
                        second_species_second_half = specimen_1[2][cros_position[0]:]

                        specimen_0[2] = np.concatenate((second_species_first_half, first_species_second_half))
                        specimen_1[2] = np.concatenate((first_species_first_half, second_species_second_half))

                    self.pop['mate_cros'].append(specimen_0)
                    self.pop['mate_cros'].append(specimen_1)

        return 1

    def ga_fitnes_function(self):

        if self.pop['num_iteration'] == 0:
            keys_list = ['pop']
        else:
            keys_list = ['pop', 'mate_mut', 'mate_cros']
        for key in keys_list:
            for i in tqdm(range(len(self.pop[key]))):
                species = self.pop[key][i][2]
                self.ep.get_data(species)

                self.add_data_to_pop(self.ep, key, i)


                fitnes_fun = self.ga_parameters['penalty_function_weights'][0] * self.error_penalty_function(self.ep) +\
                    self.ga_parameters['penalty_function_weights'][1] * self.head_penalty_function(self.ep) +\
                    self.ga_parameters['penalty_function_weights'][2] * self.flow_penalty_function(self.ep) +\
                    self.ga_parameters['penalty_function_weights'][3] * self.tank_penalty_function(self.ep) +\
                    self.ga_parameters['penalty_function_weights'][4] * self.tank_final_state_penalty_function(self.ep)

                self.pop[key][i][0] = fitnes_fun
        return self.pop

    def ga_selection(self):

        if self.ga_parameters['parent_selection_type'] == 'first_n':
            keys_list = ['pop', 'mate_mut', 'mate_cros']
            tuple_list = []
            for key in keys_list:
                for i in range(len(self.pop[key])):
                    tuple_list.append(tuple((self.pop[key][i][0],self.pop[key][i][1],self.pop[key][i][2],\
                                             self.pop[key][i][3],self.pop[key][i][4],self.pop[key][i][5],\
                                             self.pop[key][i][6])))

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


    def add_data_to_pop(self,epa,population_name,specimen_number):

        self.pop[population_name][specimen_number].append(epa.data['error_output'][0])

        temp_head = []
        for i in range(len(epa.sw_parameters['mes_nodes_names']) - len(epa.sw_parameters['tanks_names'])):
            temp_head.append(epa.data['head_output_' + epa.sw_parameters['mes_nodes_names'][i]][0])

        self.pop[population_name][specimen_number].append(temp_head)

        temp_tank = []
        for i in range(len(epa.sw_parameters['mes_nodes_names']) - len(epa.sw_parameters['tanks_names']),
                       len(epa.sw_parameters['mes_nodes_names'])):
            temp_tank.append(epa.data['head_output_' + epa.sw_parameters['mes_nodes_names'][i]][0])

        self.pop[population_name][specimen_number].append(temp_tank)

        temp_flow = []
        for i in range(len(epa.sw_parameters['mes_links_names'])):
            temp_flow.append(epa.data['flow_output_' + epa.sw_parameters['mes_links_names'][i]][0])
        self.pop[population_name][specimen_number].append(temp_flow)

    def ga_stop_criterion(self):

        output = 1

        if self.pop['num_iteration'] > self.ga_parameters['stop_criterion_min_generations']:

            diff_fitness_function = []

            last_solutions = self.pop['best_solutions'][-self.ga_parameters['stop_criterion_min_generations']:]

            for i in range(self.ga_parameters['stop_criterion_min_generations']-1, 0, -1):

                diff_fitness_function.append(abs(last_solutions[i - 1] - last_solutions[i]))

            if sum(diff_fitness_function) < self.ga_parameters['stop_criterion_min_change'] and \
                sum(self.ep.data['error_output'][0]) < 0:

                output = 0

        return output

def save_variabele_space_to_file(file_name, input_dictionary):


    key = 'last_best_solution'
    temp_date = []
    temp_date.append(f"{key}_function_value = {input_dictionary['last_best_solution'][0]}\n")
    temp_date.append(f"{key}_integer_gene = {str(input_dictionary['last_best_solution'][1]).replace(' ',';')}\n")

    temp_float_gene = str(list(input_dictionary['last_best_solution'][2])).replace(',',';')
    temp_date.append(f"{key}_float_gene = {temp_float_gene}\n")


    with open(file_name,'w') as f:
        f.writelines(temp_date)
        f.close()

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
            'penalty_function_weights': [1e9, 100, 100, 1000, 100],
            'num_generations': 200,
            'num_specimen': 150,
            'num_float_genes': sw_parameters['num_pumps'] * sw_parameters['time_duration_h'],
            'pop_int_range_low': 0,
            'pop_int_range_high': 1,
            'num_int_genes': 0,
            'parent_selection_type': 'first_n',  # rws, sus, random, tournament
            'crossover_percent_probability': 40,
            'crossover_type': 'single_point',  # two_points , uniform , scattered
            'mutation_type': 'random',  # random, inversion , scramble ,
            'specialize_operators_error_delta_value': 0.25,
            'specialize_operators_head_delta_value': 0.2,
            'mutation_percent_genes': 20,
            'mutation_percent_probability': 90,
            'pop_float_range_low': 0.5,
            'pop_float_range_high': 1,
            'stop_criterion_min_generations': 10,
            'stop_criterion_min_change': 1e-2,
        }

        ga = GA_function(ga_parameters , sw_parameters)

        ga.ga_initialization()

        ga.ga_fitnes_function()

        for i in range(ga.ga_parameters['num_generations']):

            ga.ga_mutation()

            ga.ga_crossover()

            ga.ga_fitnes_function()

            ga.ga_selection()

            print(f"Iteration number: {ga.pop['num_iteration']}, Mean population value: {ga.pop['mean_pop_value'][-1]},"
                  f" Best solution: {ga.pop['last_best_solution'][0]}")

            save_variabele_space_to_file('last_best_pop.txt', ga.pop)

            if ga.ga_stop_criterion() == 0:

                break

