import copy
import random
from statistics import mean

import numpy as np

from SW_Parameters import sw_par
from GA_Parameters import ga_par

import time

from  multiprocessing import Pool

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Get_data_from_DWDS import Epa

class GA_function:

    def __init__(self):
        self.sw_parameters = sw_par
        self.ga_parameters = ga_par
        self.ep = Epa(sw_par)
        self.pop = []
        self.plt = plt

    def ga_initialization(self):

        low_int_val = self.ga_parameters.pop['int_range_low']
        high_int_val = self.ga_parameters.pop['int_range_high']
        low_float_val = self.ga_parameters.pop['float_range_low']
        high_float_val = self.ga_parameters.pop['float_range_high']

        self.pop = {'pop': []}
        self.pop['num_iteration'] = 0
        self.pop['best_solutions'] = []
        self.pop['mean_pop_value'] = []
        self.pop['last_best_solution'] = []

        for i in range(self.ga_parameters.number['specimen']):
            specimen = [[] for i in range(8)]

            if self.ga_parameters.number['int_genes'] > 0:
                temp_int_specimen = []
                for j in range(len(self.ga_parameters.number['int_genes'])):
                    temp_int_specimen.append(random.randint(low_int_val, high_int_val + 1))

                specimen[1] = temp_int_specimen

            if self.ga_parameters.number['float_genes'] > 0:
                temp_float_specimen = []
                for j in range(self.ga_parameters.number['float_genes']):
                    temp_float_specimen.append((high_float_val - low_float_val) * random.random() + low_float_val)
                specimen[2] = temp_float_specimen
            self.pop['pop'].append(specimen)

    def penalty_function_tank_final_state(self, epa):

        for i in range(len(epa.sw_parameters.mes['nodes_names']) - len(epa.sw_parameters.tanks['names']),
                       len(epa.sw_parameters.mes['nodes_names'])):
            tank_inital_state = epa.data['head_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][0]
            tank_final_state = epa.data['head_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][-1]

        return abs(tank_inital_state - tank_final_state)

    def penalty_function_head_limits(self, epa):

        head_penalty = 0

        for j in range(epa.sw_parameters.time['duration_h']):
            for i in range(len(epa.sw_parameters.mes['nodes_names']) - len(epa.sw_parameters.tanks['names'])):
                cross_values = -epa.data['head_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][j] + \
                               epa.sw_parameters.level['min_head'][i]
                if cross_values > 0:
                    head_penalty += abs(cross_values)

                cross_values = epa.data['head_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][j] - \
                               epa.sw_parameters.level['max_head'][i]
                if cross_values > 0:
                    head_penalty += abs(cross_values)

        return head_penalty

    def penalty_function_acceletarion_pump_head(self, epa):

        acceleration_head_penalty = 0
        max_delta_pump_head = self.ga_parameters.specialize_operators['max_delta_head']

        for t in range(epa.sw_parameters.time['duration_h']-1):
            for i in self.sw_parameters.mes_head_index['pump']:
                delta_head = abs(epa.data['head_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][t] - \
                               epa.data['head_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][t+1])

                if delta_head > max_delta_pump_head:
                    acceleration_head_penalty += delta_head
        return acceleration_head_penalty

    def penalty_function_tank_limits(self, epa):
        tank_penalty = 0

        for j in range(epa.sw_parameters.time['duration_h']):
            for i in range(len(epa.sw_parameters.tanks['names'])):
                cross_values = -epa.data['tank_output_' + epa.sw_parameters.tanks['names'][i]][0][j] + \
                               epa.sw_parameters.level['min_tank'][i]
                if cross_values > 0:
                    tank_penalty += abs(cross_values)

                cross_values = epa.data['tank_output_' + epa.sw_parameters.tanks['names'][i]][0][j] - \
                               epa.sw_parameters.level['max_tank'][i]
                if cross_values > 0:
                    tank_penalty += abs(cross_values)

        return tank_penalty

    def penalty_function_flow_limits(self, epa):

        flow_penalty = 0
        for j in range(epa.sw_parameters.time['duration_h']):
            for i in range(len(epa.sw_parameters.mes['links_names'])):
                cross_values = -epa.data['flow_output_' + epa.sw_parameters.mes['links_names'][i]][0][j] + \
                               epa.sw_parameters.level['min_flow'][i]

                if cross_values > 0:
                    flow_penalty += abs(cross_values)

            for i in range(len(epa.sw_parameters.mes['links_names'])):
                cross_values = epa.data['flow_output_' + epa.sw_parameters.mes['links_names'][i]][0][j] - \
                               epa.sw_parameters.level['max_flow'][i]
                if cross_values > 0:
                    flow_penalty += abs(cross_values)
        return flow_penalty

    def penalty_function_error(self,epa):

        error_values = epa.data['error_output'][0]
        return sum(error_values)

    def penalty_function_energy(self,epa):

        energy_penalty = 0

        for name in epa.sw_parameters.pumps['names']:
            for t in range(self.sw_parameters.time['duration_h']):
                energy_penalty += self.sw_parameters.energy_taryf[t] * \
                                  epa.data['energy_output_'+name][0][t]

        return energy_penalty

    def ga_mutation(self):

        self.pop['mate_mut'] = []

        self.specimen = copy.deepcopy(self.pop['pop'])

        if random.random() < 0.5:
            if self.ga_parameters.mutation['type'] == 'random':
                self.ga_mutation_random()
        else:
            #if sum(self.ep.data['error_output'][0]) > 0:
            self.ga_mutiation_SGO_error()
            #else:
            self.ga_mutiation_SGO_head()
            self.ga_mutiation_SGO_flow()
            self.ga_mutiation_SGO_end_tank_level()
            self.ga_mutiation_SGO_acceleration_pump_head()

    def ga_mutation_random(self):

        low_int_val = self.ga_parameters.pop['int_range_low']
        high_int_val = self.ga_parameters.pop['int_range_high']
        low_float_val = self.ga_parameters.pop['float_range_low']
        high_float_val = self.ga_parameters.pop['float_range_high']


        for i in range(len(self.specimen)):
            if random.random() < self.ga_parameters.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if self.ga_parameters.number['int_genes'] > 0:

                    num_genes = int(self.ga_parameters.number['int_genes']*\
                                    self.ga_parameters.mutation['percent_genes']/100)

                    genes_position =[]
                    for j in range(num_genes):
                        genes_position.append(random.randint(0,self.ga_parameters.number['int_genes']-1))

                    for ge_pos in genes_position:
                        species[1][ge_pos] = random.randint(low_int_val, high_int_val+1)

                    self.pop['mate_mut'].append(species)

                if self.ga_parameters.number['float_genes'] > 0:
                    if random.random() < self.ga_parameters.mutation['percent_probability']/100:

                       num_genes = int(
                           self.ga_parameters.number['float_genes'] * self.ga_parameters.mutation['percent_genes'] / 100)

                       genes_position = []
                       for j in range(num_genes):
                           genes_position.append(random.randint(0, self.ga_parameters.number['float_genes']-1))

                       for ge_pos in genes_position:
                           if random.random()<0.5:
                                species[2][ge_pos] += random.random() * self.ga_parameters.mutation['random_delta']
                           else:
                                species[2][ge_pos] -= random.random() * self.ga_parameters.mutation['random_delta']

                           if species[2][ge_pos] > high_float_val:
                               species[2][ge_pos] = high_float_val
                           if species[2][ge_pos] < low_float_val:
                               species[2][ge_pos] = low_float_val
                    self.pop['mate_mut'].append(species)

    def ga_mutiation_SGO_error(self):

        time_dur = self.sw_parameters.time['duration_h']
        for i in range(len(self.specimen)):
            if random.random() < self.ga_parameters.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if self.ga_parameters.number['float_genes'] > 0:
                    for j in range(len(species[3])):
                        error_value = species[3][j]
                        flow_pump_karolewo = species[6][0][j]
                        flow_pump_funka = species[6][1][j]
                        flow_pump_plac = species[6][2][j]
                        if error_value > 4 and flow_pump_karolewo == 0 and flow_pump_plac > 0:
                            species[2][j + time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                'error_delta_value']
                            species[2][j + 2 * time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                'error_delta_value']

                        if error_value > 4 and flow_pump_funka == 0:
                            if random.random() < 0.5:
                                species[2][j] += random.random() * self.ga_parameters.specialize_operators[
                                    'error_delta_value']*2
                            else:
                                species[2][j] = 1


                        if error_value > 4 and flow_pump_plac == 0 and flow_pump_karolewo > 0:
                            species[2][j + 2 * time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                'error_delta_value']
                            species[2][j + time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                'error_delta_value']

                        if error_value > 4 and flow_pump_plac == 0 and flow_pump_karolewo == 0:
                            species[2][j + 2 * time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                'error_delta_value']
                            species[2][j + time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                'error_delta_value']

                self.pop['mate_mut'].append(self.set_specimen_bound(species))

    def ga_mutiation_SGO_head(self):

        time_dur = self.sw_parameters.time['duration_h']
        for i in range(len(self.specimen)):
            if random.random() < self.ga_parameters.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if self.ga_parameters.number['float_genes'] > 0:
                    for j in range(time_dur):
                        kk=0
                        for k in range(len(self.sw_parameters.mes['nodes_names'])):

                            head_min_lev = self.sw_parameters.level['min_head'][k]
                            head_max_level = self.sw_parameters.level['max_head'][k]

                            if k < len(self.sw_parameters.mes['nodes_names']) - len(self.sw_parameters.tanks['names']):
                                head_WMO = species[4][k][j]
                                delta_min_head = - head_min_lev + head_WMO
                                delta_max_head = - head_WMO + head_max_level
                            else:
                                tank_level = species[5][kk][j]
                                delta_min_head = - head_min_lev + tank_level
                                delta_max_head = - tank_level + head_max_level
                                kk = kk + 1

                            if k == 0:
                                if delta_min_head < 0:
                                    species[2][j] += random.random() * self.ga_parameters.specialize_operators[
                                    'head_delta_value']

                                if delta_max_head < 0:
                                    species[2][j] -= random.random() * self.ga_parameters.specialize_operators[
                                        'head_delta_value']
                            if k == 1:
                                if  delta_min_head < 0:
                                    species[2][j + time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                    'head_delta_value']
                                if delta_max_head < 0:
                                    species[2][j + time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                        'head_delta_value']

                            if k == 2:
                                if  delta_min_head < 0:
                                    species[2][j + 2 * time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                    'head_delta_value']

                                if delta_max_head < 0:
                                    species[2][j + 2 * time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                        'head_delta_value']

                            if k == 3:
                                if  delta_min_head < 0:
                                    species[2][j + time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                    'head_delta_value']
                                    species[2][j + 2 * time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                    'head_delta_value']

                                if delta_max_head < 0:
                                    species[2][j + time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                        'head_delta_value']
                                    species[2][j + 2 * time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                        'head_delta_value']

                            if k == 4:
                                if delta_min_head < 0:
                                    species[2][j + time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                    'head_delta_value']
                                    species[2][j + 2 * time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                    'head_delta_value']

                                if delta_max_head < 0:
                                    species[2][j + time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                        'head_delta_value']
                                    species[2][j + 2 * time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                        'head_delta_value']
                            #Zbiorniki
                            if k == 5:
                                if delta_min_head < 0:
                                    species[2][j] += random.random() * \
                                                     self.ga_parameters.specialize_operators['head_delta_value']
                                    species[2][j + time_dur] -= random.random() \
                                                                * self.ga_parameters.specialize_operators[
                                                                'head_delta_value']
                                    species[2][j + 2*time_dur] += random.random() \
                                                                * self.ga_parameters.specialize_operators[
                                                                'head_delta_value']
                                if delta_max_head < 0:
                                    species[2][j] -= random.random() * \
                                                     self.ga_parameters.specialize_operators['head_delta_value']
                                    species[2][j + time_dur] += random.random() \
                                                                * self.ga_parameters.specialize_operators[
                                                                    'head_delta_value']
                                    species[2][j + 2 * time_dur] -= random.random() \
                                                                    * self.ga_parameters.specialize_operators[
                                                                        'head_delta_value']

                self.pop['mate_mut'].append(self.set_specimen_bound(species))

    def ga_mutiation_SGO_acceleration_pump_head(self):

        time_dur = self.sw_parameters.time['duration_h']
        for i in range(len(self.specimen)):
            if random.random() < self.ga_parameters.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if self.ga_parameters.number['float_genes'] > 0:
                    for j in range(time_dur-1):
                        kk = 0
                        for k in self.sw_parameters.mes_head_index['pump']:

                            max_delta_pump_head = self.ga_parameters.specialize_operators['max_delta_head']

                            delta_pump_head = species[4][k][j] - species[4][k][j+1]

                            if k == 0:

                                if delta_pump_head < - max_delta_pump_head:
                                    species[2][j] += random.random() * self.ga_parameters.specialize_operators[
                                                                        'head_delta_value']
                                    species[2][j+1] -= random.random() * self.ga_parameters.specialize_operators[
                                                                        'head_delta_value']
                                if delta_pump_head > max_delta_pump_head:
                                    species[2][j] -= random.random() * self.ga_parameters.specialize_operators[
                                                                        'head_delta_value']
                                    species[2][j] += random.random() * self.ga_parameters.specialize_operators[
                                                                        'head_delta_value']
                            if k == 1:

                                if delta_pump_head < - max_delta_pump_head:
                                    species[2][j + time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                                                        'head_delta_value']
                                    species[2][j+1+ time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                                                        'head_delta_value']
                                if delta_pump_head > max_delta_pump_head:
                                    species[2][j+ time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                                                        'head_delta_value']
                                    species[2][j+ time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                                                        'head_delta_value']
                            if k == 2:

                                if delta_pump_head < - max_delta_pump_head:
                                    species[2][j + 2*time_dur] += random.random() * \
                                                                self.ga_parameters.specialize_operators[
                                                                    'head_delta_value']
                                    species[2][j + 1 + 2*time_dur] -= random.random() * \
                                                                self.ga_parameters.specialize_operators[
                                                                    'head_delta_value']
                                if delta_pump_head > max_delta_pump_head:

                                    species[2][j + 2*time_dur] -= random.random() * \
                                                                self.ga_parameters.specialize_operators[
                                                                    'head_delta_value']
                                    species[2][j + 2*time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                    'head_delta_value']

                self.pop['mate_mut'].append(self.set_specimen_bound(species))

    def ga_mutiation_SGO_end_tank_level(self):

        time_dur = self.sw_parameters.time['duration_h']
        for i in range(len(self.specimen)):
            if random.random() < self.ga_parameters.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if self.ga_parameters.number['float_genes'] > 0:

                    for j in range(len(self.sw_parameters.tanks['names'])):

                        init_level = species[5][j][0]
                        end_level = species[5][j][-1]
                        if end_level - init_level > self.ga_parameters.specialize_operators['delta_initial_end_tank_level']:
                            for k in range(time_dur - self.ga_parameters.specialize_operators['end_tank_level_horizon'],
                                           time_dur):
                                species[2][j + k] -= random.random() * self.ga_parameters.specialize_operators[
                                    'head_delta_value']

                                species[2][j + k + time_dur] += random.random() * \
                                                                self.ga_parameters.specialize_operators['head_delta_value']

                                species[2][j + k + 2 * time_dur] -= random.random() * \
                                                                    self.ga_parameters.specialize_operators['head_delta_value']

                        if end_level - init_level < - self.ga_parameters.specialize_operators['delta_initial_end_tank_level']:
                            for k in range(time_dur - self.ga_parameters.specialize_operators['end_tank_level_horizon'],
                                           time_dur):
                                species[2][j + k] += random.random() * self.ga_parameters.specialize_operators[
                                    'head_delta_value']

                                species[2][j + k + time_dur] -= random.random() * \
                                                                self.ga_parameters.specialize_operators['head_delta_value']

                                species[2][j + k + 2 * time_dur] += random.random() * \
                                                                    self.ga_parameters.specialize_operators['head_delta_value']

                self.pop['mate_mut'].append(self.set_specimen_bound(species))


    def ga_mutiation_SGO_flow(self):

        time_dur = self.sw_parameters.time['duration_h']
        for i in range(len(self.specimen)):
            if random.random() < self.ga_parameters.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if self.ga_parameters.number['float_genes'] > 0:
                    for j in range(time_dur):
                        for k in range(len(self.sw_parameters.mes['links_names'])):

                            flow_mes = species[4][k][j]

                            flow_min_level = self.sw_parameters.level['min_flow'][k]
                            flow_max_level = self.sw_parameters.level['max_flow'][k]

                            delta_min_flow = - flow_min_level + flow_mes
                            delta_max_flow = - flow_mes + flow_max_level
                            if k == 0:
                                if delta_min_flow < 0:
                                    species[2][j] += random.random() * self.ga_parameters.specialize_operators[
                                    'flow_delta_value']

                                if delta_max_flow < 0:
                                    species[2][j] -= random.random() * self.ga_parameters.specialize_operators[
                                        'flow_delta_value']
                            if k == 1:
                                if delta_min_flow < 0:
                                    species[2][j + time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                    'flow_delta_value']
                                if delta_max_flow < 0:
                                    species[2][j + time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                        'flow_delta_value']

                            if k == 2:
                                if delta_min_flow < 0:
                                    species[2][j + 2 * time_dur] += random.random() * self.ga_parameters.specialize_operators[
                                    'flow_delta_value']

                                if delta_max_flow < 0:
                                    species[2][j + 2 * time_dur] -= random.random() * self.ga_parameters.specialize_operators[
                                        'flow_delta_value']

                self.pop['mate_mut'].append(self.set_specimen_bound(species))

    def set_specimen_bound(self, species):
        low_float_val = self.ga_parameters.pop['float_range_low']
        high_float_val = self.ga_parameters.pop['float_range_high']

        for t in range(len(species)):
            if species[2][t] > high_float_val:
                species[2][t] = high_float_val
            if species[2][t] < low_float_val:
                species[2][t] = low_float_val

        return species

    def ga_crossover(self):

        #low_int_val = self.ga_parameters['pop_int_range_low']
        #high_int_val = self.ga_parameters['pop_int_range_high']
        #low_float_val = self.ga_parameters['pop_float_range_low']
        #high_float_val = self.ga_parameters['pop_float_range_high']

        temp_cross = copy.deepcopy(self.pop['pop'])

        for te in self.pop['mate_mut']:
            temp_cross.append(te)

        self.pop['mate_cros'] = []
        if self.ga_parameters.crossover['type'] == 'single_point':
            for i in range(self.ga_parameters.number['specimen']):
                if random.random() < self.ga_parameters.crossover['percent_probability'] / 100:

                    specimen_to_mutation = [random.randint(1, len(temp_cross)-1)]
                    specimen_to_mutation.append(random.randint(1, len(temp_cross)-1))

                    specimen_0 = list(copy.deepcopy(temp_cross[specimen_to_mutation[0]]))
                    specimen_1 = list(copy.deepcopy(temp_cross[specimen_to_mutation[1]]))

                    if self.ga_parameters.number['int_genes'] > 0:
                        cros_position = random.randint(1, self.ga_parameters.number['int_genes'], size=(1,))

                        first_species_first_half = specimen_0[1][:cros_position[0]]
                        first_species_second_half = specimen_0[1][cros_position[0]:]
                        second_species_first_half = specimen_1[1][:cros_position[0]]
                        second_species_second_half = specimen_1[1][cros_position[0]:]

                        specimen_0[1] = second_species_first_half + first_species_second_half
                        specimen_1[1] = first_species_first_half + second_species_second_half

                    if self.ga_parameters.number['float_genes'] > 0:
                        cros_position = [random.randint(1, self.ga_parameters.number['float_genes'])]

                        first_species_first_half = specimen_0[2][:cros_position[0]]
                        first_species_second_half = specimen_0[2][cros_position[0]:]
                        second_species_first_half = specimen_1[2][:cros_position[0]]
                        second_species_second_half = specimen_1[2][cros_position[0]:]

                        specimen_0[2] = second_species_first_half + first_species_second_half
                        specimen_1[2] = first_species_first_half + second_species_second_half

                    self.pop['mate_cros'].append(specimen_0)
                    self.pop['mate_cros'].append(specimen_1)

        return 1

    def add_data_to_pop(self,epa,population_name,specimen_number):

        #3
        self.pop[population_name][specimen_number][3] = epa.data['error_output'][0]

        temp_head = []#4
        for i in range(len(epa.sw_parameters.mes['nodes_names']) - len(epa.sw_parameters.tanks['names'])):
            temp_head.append(epa.data['head_output_' + epa.sw_parameters.mes['nodes_names'][i]][0])
        self.pop[population_name][specimen_number][4] = temp_head

        temp_tank = []#5
        for i in range(len(epa.sw_parameters.tanks['names'])):
            temp_tank.append(epa.data['tank_output_' + epa.sw_parameters.tanks['names'][i]][0])

        self.pop[population_name][specimen_number][5] = temp_tank

        temp_flow = []#6
        for ln in epa.sw_parameters.mes['links_names']:
            temp_flow.append(epa.data['flow_output_' + ln][0])
        self.pop[population_name][specimen_number][6] = temp_flow

    def ga_fitnes_function2(self):

        if self.pop['num_iteration'] == 0:
            keys_list = ['pop']
        else:
            keys_list = ['pop', 'mate_mut', 'mate_cros']

        pomoc = []


        for k in range(len(keys_list)):
            with Pool() as pool:
                items = [('pop', 1), ('pop',2), ('pop',3)]
                print(items)
                for res in pool.starmap(prarallel_com, items):
                    print(res)


        self.pop['num_iteration'] += 1
        return self.pop

    def ga_fitnes_function(self):

        if self.pop['num_iteration'] == 0:
            keys_list = ['pop']
        else:
            keys_list = ['pop', 'mate_mut', 'mate_cros']

        pomoc =[]
        for k in range(len(keys_list)):
            for i in (range(len(self.pop[keys_list[k]]))):
                species = self.pop[keys_list[k]][i][2]

                self.ep.get_data(species)

                self.add_data_to_pop(self.ep, keys_list[k], i)

                penalties = {}
                penalties['error'] = self.ga_parameters.penalty_function_weights[0][0] * \
                                     self.penalty_function_error(self.ep)
                penalties['energy'] = self.ga_parameters.penalty_function_weights[0][1] * \
                                      self.penalty_function_energy(self.ep)
                penalties['head'] = self.ga_parameters.penalty_function_weights[0][2] * \
                                    self.penalty_function_head_limits(self.ep)
                penalties['flow'] = self.ga_parameters.penalty_function_weights[0][3] * \
                                    self.penalty_function_flow_limits(self.ep)
                penalties['tank'] = self.ga_parameters.penalty_function_weights[0][4] * \
                                    self.penalty_function_tank_limits(self.ep)
                penalties['tank_final'] = self.ga_parameters.penalty_function_weights[0][5] * \
                                          self.penalty_function_tank_final_state(self.ep)
                penalties['pump_acc'] = self.ga_parameters.penalty_function_weights[0][6] * \
                                        self.penalty_function_acceletarion_pump_head(self.ep)

                fitnes_fun =  penalties['error'] +  penalties['energy'] + penalties['head'] + penalties['flow'] + \
                              penalties['tank'] + penalties['tank_final'] + penalties['pump_acc']

                self.pop[keys_list[k]][i][0] = fitnes_fun
                self.pop[keys_list[k]][i][7] = penalties

        self.pop['num_iteration'] += 1

        return self.pop

    def ga_selection(self):
        if self.ga_parameters.selection[0]['type'] == 'first_n':
            keys_list = ['pop', 'mate_mut', 'mate_cros']
            ii=0
            pop_list=[]
            pomoc=[]
            for key in keys_list:
                for i in range(len(self.pop[key])):
                    ii+=1
                    pomoc.append([self.pop[key][i][0], key])
                    pop_list.append(self.pop[key][i])
            self.pop['pop'] = sorted(pop_list)[:int(self.ga_parameters.number['specimen']*\
                                                    self.ga_parameters.selection[0]['n_type_precent'])]



        '''if self.ga_parameters['parent_selection_type'] == 'roulete':
            keys_list = ['pop', 'mate_mut', 'mate_cros']
            ii=0
            pop_list=[]
            pomoc=[]
            for key in keys_list:
                for i in range(len(self.pop[key])):'''


        return 0

    def ga_statistics(self):

        self.pop['last_best_solution'] = self.pop['pop'][0]

        self.pop['best_solutions'].append(self.pop['last_best_solution'][0])

        population_fun_fitnes = []

        for i in range(ga_par.number['specimen']):
            population_fun_fitnes.append(self.pop['pop'][i][0])

        self.pop['mean_pop_value'].append(mean(population_fun_fitnes))

    def ga_stop_criterion(self):

        output = 1

        if self.pop['num_iteration'] > self.ga_parameters.stop_criterion['min_generations']:

            diff_fitness_function = []

            last_solutions = self.pop['best_solutions'][-self.ga_parameters.stop_criterion['min_generations']:]

            for i in range(self.ga_parameters.stop_criterion['min_generations']-1, 0, -1):
                diff_fitness_function.append(abs(last_solutions[i - 1] - last_solutions[i]))

            if sum(diff_fitness_function) < self.ga_parameters.stop_criterion['min_change'] and \
                sum(self.ep.data['error_output'][0]) < 0:

                output = 0

        return output

    def print_statistics(self):

        self.ga_statistics()
        self.print_name = f"It. num.: {ga.pop['num_iteration'] - 1}, "+\
                          f"Mean: {ga.pop['mean_pop_value'][-1]:.0f},"+ \
                          f"Best: {ga.pop['last_best_solution'][0]:.0f}, "+\
                          f"Energy: {ga.pop['last_best_solution'][7]['energy']:.0f}, "+\
                          f"Head: {ga.pop['last_best_solution'][7]['head']:.0f}, "+\
                          f"Flow: {ga.pop['last_best_solution'][7]['flow']:.0f}, "+\
                          f"Tank: {ga.pop['last_best_solution'][7]['tank']:.0f}, "+\
                          f"Tank final: {ga.pop['last_best_solution'][7]['tank_final']:.0f}, " + \
                          f"Head acc: {ga.pop['last_best_solution'][7]['pump_acc']:.0f}"
        print('--------------------------------------------------------------------------------------------')
        print(self.print_name)

        #print(f"error: {ga.pop['last_best_solution'][3]}")
        #print(f"head 1: {ga.pop['last_best_solution'][4][0]}")
        #print(f"head 2: {ga.pop['last_best_solution'][4][1]}")
        #print(f"head 55: {ga.pop['last_best_solution'][4][2]}")
        #print(f"head 83: {ga.pop['last_best_solution'][4][3]}")
        #print(f"head 158: {ga.pop['last_best_solution'][4][4]}")
        #print(f"tank: {ga.pop['last_best_solution'][5][0]}")
        #print(f"flow F: {ga.pop['last_best_solution'][6][0]}")
        #print(f"flow K: {ga.pop['last_best_solution'][6][1]}")
        #print(f"flow P: {ga.pop['last_best_solution'][6][2]}")

    def plot_mes(self, axis, type, which, location):

        x = np.arange(0, self.sw_parameters.time['duration_h'])

        if type == 4:
            mes_name = 'nodes'
            name = 'head'
            ylabel = 'Head [m]'
        if type ==5:
            mes_name = 'tank'
            name = 'tank'
            ylabel = 'Head [m]'
        if type == 6:
            mes_name = 'links'
            name = 'flow'
            ylabel = 'Flow [m^3/h]'

        upper_bound = [self.sw_parameters.level['max_'+name][which[0]]]*self.sw_parameters.time['duration_h']
        lower_bound = [self.sw_parameters.level['min_'+name][which[0]]]*self.sw_parameters.time['duration_h']
        self.ax[axis[0], axis[1]].clear()
        color = ['b','m','g']
        for i in range(len(which)):
            self.ax[axis[0], axis[1]].plot(x, self.pop['last_best_solution'][type][which[i]],
                                            color=color[i], linewidth=1)
        self.ax[axis[0], axis[1]].plot(x, lower_bound, color='r', linewidth=1)
        self.ax[axis[0], axis[1]].plot(x, upper_bound, color='r', linewidth=1)
        self.ax[axis[0], axis[1]].set_ylabel(ylabel, fontsize=9)
        self.ax[axis[0], axis[1]].set_xlabel("Time [h]", fontsize=9)


        labels =[]
        names =''
        for i in range(len(which)):
            labels.append(self.sw_parameters.mes[mes_name + '_names'][which[i]])
            if i < len(which)-1:
                names += self.sw_parameters.mes[mes_name +'_names'][which[i]]+', '
            else:
                names += self.sw_parameters.mes[mes_name +'_names'][which[i]]

        self.ax[axis[0], axis[1]].set_title(name + ' at '+ names,
                                            fontsize = 10)
        self.ax[axis[0], axis[1]].xaxis.set_tick_params(labelsize=8)
        self.ax[axis[0], axis[1]].yaxis.set_tick_params(labelsize=6)

        self.ax[axis[0], axis[1]].legend(loc = location, fontsize='8',
                                         labels = labels)

    def ga_init_plot(self):
        self.plt.ion()
        self.plot = {}
        self.plot['best_sol'] = self.pop['best_solutions']
        self.plot['mean_sol'] = [self.pop['mean_pop_value']]
        self.plot['x'] = range(len(self.plot['best_sol']))
        figure, ax = plt.subplots(nrows=2, ncols=4, figsize = (14,10))
        mngr = plt.get_current_fig_manager()
        mngr.window.geometry("+100+0")
        self.figure = figure
        self.ax = ax

    def plot_fitnes_values(self,place,i):
        x =np.arange(0, i+1)

        self.ax[place[0], place[1]].plot(x, self.plot['best_sol'], label = 'Best specimen',
                                         color = 'r', linewidth = 1)
        self.ax[place[0], place[1]].plot(x, self.plot['mean_sol'][0], label = 'Mean of specimens',
                                         color = 'g', linewidth = 1)
        self.ax[place[0], place[1]].set_title("Mean and best solution", fontsize='8')
        self.ax[place[0], place[1]].set_ylabel("Mean and best", fontsize='8')
        self.ax[place[0], place[1]].set_yscale('log')
        self.ax[place[0], place[1]].set_xlabel("Iteration", fontsize='8')
        self.ax[place[0], place[1]].xaxis.set_tick_params(labelsize=8)
        self.ax[place[0], place[1]].yaxis.set_tick_params(labelsize=6)
        self.ax[place[0], place[1]].legend(loc="upper right", fontsize='8',
                                         labels = ['Best value', 'Mean_value'])

    def print_values_on_plot(self, epa):
        #Rectangle((20, 0), 50, 50, facecolor='red')
        error = ''
        i = 0
        for er in epa['pop'][0][3]:
            if er == 4:
                error += str(i) + ', '
            i += 1
        self.figure.text(0.01, 0.95, ' ' * 800, fontsize=12, bbox=dict(facecolor='white', edgecolor='white'))
        self.figure.text(0.01, 0.95, self.print_name, fontsize=12, bbox ={'facecolor':'white'})
        self.figure.text(0.01, 0.92, ' '*100, fontsize=12, bbox=dict(facecolor= 'white', edgecolor = 'white'))
        self.figure.text(0.01, 0.92, error, fontsize=12, bbox={'facecolor': 'white'})
    def ga_plot(self,i, epa):

        self.plot_fitnes_values([0,0],i)
        self.plot_mes([0, 1], 4, [0],'upper right')
        self.plot_mes([0, 2], 4, [1, 2], 'upper right')
        self.plot_mes([0, 3], 4, [3], 'lower right')
        self.plot_mes([1, 0], 4, [4], 'lower right')
        self.plot_mes([1, 1], 5, [0], 'lower right')
        self.plot_mes([1, 2], 6, [0, 1], 'lower right')
        self.plot_mes([1, 3], 6, [2], 'lower right')
        self.print_values_on_plot(epa)
        self.figure.canvas.draw()

        self.figure.canvas.flush_events()

        time.sleep(0.1)

def prarallel_com(key_list, i):
    print(input)

    species = self.pop[key_list][i][2]

    self.ep.get_data(species)

    self.add_data_to_pop(self.ep, key_list, i)

    penalties = {}
    penalties['error'] = self.ga_parameters.penalty_function_weights[0][
                             0] * self.error_penalty_function(self.ep)
    penalties['energy'] = self.ga_parameters.penalty_function_weights[0][
                              1] * self.energy_penalty_function(self.ep)
    penalties['head'] = self.ga_parameters.penalty_function_weights[0][2] * self.head_penalty_function(
        self.ep)
    penalties['flow'] = self.ga_parameters.penalty_function_weights[0][3] * self.flow_penalty_function(
        self.ep)
    penalties['tank'] = self.ga_parameters.penalty_function_weights[0][4] * self.tank_penalty_function(
        self.ep)
    penalties['tank_final'] = self.ga_parameters.penalty_function_weights[0][
                                  5] * self.tank_final_state_penalty_function(self.ep)

    fitnes_fun = penalties['error'] + penalties['energy'] + penalties['head'] + penalties['flow'] + \
                 penalties['tank'] + penalties['tank_final']

    self.pop[key_list][i][0] = fitnes_fun
    self.pop[key_list][i][7] = penalties

def save_variabele_space_to_file(file_name, input_dictionary):

    key = 'last_best_solution'
    temp_date = []
    temp_date.append(f"{key}_function_value = {input_dictionary[key][0]}\n")
    temp_date.append(f"{key}_integer_gene = {str(input_dictionary[key][1]).replace(' ',';')}\n")

    temp_float_gene = str(list(input_dictionary[key][2])).replace(',',';')
    temp_date.append(f"{key}_float_gene = {temp_float_gene}\n")


    with open(file_name,'w') as f:
        f.writelines(temp_date)
        f.close()

if __name__ == '__main__':
        print('\nTo jest biblioteka pomocna przy algorytmie genetycznym')

        ga = GA_function()

        ga.ga_initialization()

        ga.ga_fitnes_function()
        #ga.ga_fitnes_function2()

        ga.ga_statistics()

        ga.ga_init_plot()

        for i in range(1,ga.ga_parameters.number['generations']+1):

            ga.ga_mutation()

            ga.ga_crossover()

            ga.ga_fitnes_function()

            ga.ga_selection()

            save_variabele_space_to_file('last_best_pop.txt', ga.pop)

            ga.print_statistics()

            ga.ga_plot(i, ga.pop)

            if ga.ga_stop_criterion() == 0:

                break