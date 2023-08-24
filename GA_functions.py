import copy
import random
from statistics import mean

import numpy as np

from SW_Parameters import sw_par
from GA_Parameters import ga_par

import time

from  multiprocessing import Pool

import matplotlib.pyplot as plt
from Get_data_from_DWDS import Epa

class GA_function:

    def __init__(self):
        self.ep = Epa(sw_par)
        self.pop = []
        self.plt = plt
        self.init_sw_parameters()

    def init_sw_parameters(self):
        self.ep.open_epanet()
        self.ep.get_link_index()
        self.ep.get_node_index()
        self.ep.get_node_pattern()
        self.ep.get_node_pattern_names()
        self.ep.get_pumps_pattern_index()
        self.ep.get_pattern_values()

        self.ep.close_epanet()

    def ga_initialization(self):

        low_int_val = ga_par.pop['int_range_low']
        high_int_val = ga_par.pop['int_range_high']
        low_float_val = ga_par.pop['float_range_low']
        high_float_val = ga_par.pop['float_range_high']

        self.pop = {'pop': []}
        self.pop['num_iteration'] = 0
        self.pop['best_solutions'] = []
        self.pop['mean_pop_value'] = []
        self.pop['last_best_solution'] = []

        for i in range(ga_par.number['specimen']):
            specimen = [[] for i in range(9)]

            if ga_par.number['int_genes'] > 0:
                temp_int_specimen = []
                for j in range(len(ga_par.number['int_genes'])):
                    temp_int_specimen.append(random.randint(low_int_val, high_int_val + 1))

                specimen[1] = temp_int_specimen

            if ga_par.number['float_genes'] > 0:
                temp_float_specimen = []
                for j in range(ga_par.number['float_genes']):
                    temp_float_specimen.append((high_float_val - low_float_val) * random.random() + low_float_val)
                specimen[2] = temp_float_specimen
            self.pop['pop'].append(specimen)

    def penalty_function_tank_final_state(self, epa):

        penalty_init_end_tank_level = 0
        for i in sw_par.mes_pressure_index['tank']:
            tank_inital_state = epa.data['pressure_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][0]
            tank_final_state = epa.data['pressure_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][-1]

            delta_init_end_tank = abs(tank_final_state - tank_inital_state)

            if delta_init_end_tank > ga_par.specialize_operators['delta_initial_end_tank_level']:
                penalty_init_end_tank_level += delta_init_end_tank

        return (10*abs(penalty_init_end_tank_level))**2

    def penalty_function_pressure_limits(self, epa):

        pressure_penalty = 0

        for j in range(epa.sw_parameters.time['duration_h']):
            for i in range(len(epa.sw_parameters.mes['nodes_names']) - len(epa.sw_parameters.tanks['names'])):
                cross_values = -epa.data['pressure_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][j] + \
                               epa.sw_parameters.level['min_pressure'][i]
                if cross_values > 0:
                    pressure_penalty += (10*abs(cross_values))**2

                cross_values = epa.data['pressure_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][j] - \
                               epa.sw_parameters.level['max_pressure'][i]
                if cross_values > 0:
                    pressure_penalty += (10*abs(cross_values))**2

        return pressure_penalty

    def penalty_function_acceletarion_pump_pressure(self, epa):

        acceleration_pressure_penalty = 0
        max_delta_pump_pressure = ga_par.specialize_operators['max_delta_pump_pressure']

        for t in range(epa.sw_parameters.time['duration_h']-1):
            for i in sw_par.mes_pressure_index['pump']:
                delta_pressure = abs(epa.data['pressure_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][t] - \
                               epa.data['pressure_output_' + epa.sw_parameters.mes['nodes_names'][i]][0][t+1])

                if delta_pressure > max_delta_pump_pressure:
                    acceleration_pressure_penalty += (10*delta_pressure)**2
        return acceleration_pressure_penalty

    def penalty_function_acceletarion_pump_speed(self, epa):

        acceleration_speed_penalty = 0
        max_delta_pump_speed = ga_par.specialize_operators['max_delta_pump_speed']

        for t in range(sw_par.time['duration_h']-1):
            for i in range(len(sw_par.mes['links_names'])):
                delta_speed = abs(epa.data['flow_output_' + epa.sw_parameters.mes['links_names'][i]][0][t] - \
                               epa.data['flow_output_' + epa.sw_parameters.mes['links_names'][i]][0][t+1])

                if delta_speed > max_delta_pump_speed:
                    acceleration_speed_penalty += (10*delta_speed)**2
        return acceleration_speed_penalty

    def penalty_function_acceletarion_pump_flow(self, epa):

        acceleration_flow_penalty = 0
        max_delta_pump_flow = ga_par.specialize_operators['max_delta_pump_flow']

        for t in range(epa.sw_parameters.time['duration_h']-1):
            for i in range(len(epa.sw_parameters.mes['links_names'])):
                delta_flow = abs(epa.data['flow_output_' + epa.sw_parameters.mes['links_names'][i]][0][t] - \
                               epa.data['flow_output_' + epa.sw_parameters.mes['links_names'][i]][0][t+1])

                if delta_flow > max_delta_pump_flow:
                    acceleration_flow_penalty += (10*delta_flow)**2
        return acceleration_flow_penalty

    def penalty_function_tank_limits(self, epa):
        tank_penalty = 0

        for j in range(epa.sw_parameters.time['duration_h']):
            for i in range(len(epa.sw_parameters.tanks['names'])):
                cross_values = -epa.data['tank_output_' + epa.sw_parameters.tanks['names'][i]][0][j] + \
                               epa.sw_parameters.level['min_tank'][i]
                if cross_values > 0:
                    tank_penalty += (10*abs(cross_values))**2

                cross_values = epa.data['tank_output_' + epa.sw_parameters.tanks['names'][i]][0][j] - \
                               epa.sw_parameters.level['max_tank'][i]
                if cross_values > 0:
                    tank_penalty += (10*abs(cross_values))**2

        return tank_penalty

    def penalty_function_flow_limits(self, epa):

        flow_penalty = 0
        for j in range(epa.sw_parameters.time['duration_h']):
            for i in range(len(epa.sw_parameters.mes['links_names'])):
                cross_values = -epa.data['flow_output_' + epa.sw_parameters.mes['links_names'][i]][0][j] + \
                               epa.sw_parameters.level['min_flow'][i]

                if cross_values > 0:
                    flow_penalty += (10*abs(cross_values))**2

            for i in range(len(epa.sw_parameters.mes['links_names'])):
                cross_values = epa.data['flow_output_' + epa.sw_parameters.mes['links_names'][i]][0][j] - \
                               epa.sw_parameters.level['max_flow'][i]
                if cross_values > 0:
                    flow_penalty += (10*abs(cross_values))**2
        return flow_penalty

    def penalty_function_error(self,epa):

        error_values = epa.data['error_output'][0]

        for j in range(epa.sw_parameters.time['duration_h']):
            k = 0
            for name in epa.sw_parameters.mes['links_names']:
                if epa.data['flow_output_' + name][0][j] <= 0 :
                    k += 1
            if error_values[j] > 0:
                error_values[j] = error_values[j] #** k

        return sum(error_values)

    def penalty_function_energy(self,epa):

        energy_penalty = 0

        for name in epa.sw_parameters.pumps['names']:
            for t in range(sw_par.time['duration_h']):
                energy_penalty += sw_par.energy_taryf[t] * \
                                  epa.data['energy_output_'+name][0][t]

        return energy_penalty

    def ga_mutation(self):

        self.pop['mate_mut'] = []

        self.specimen = copy.deepcopy(self.pop['pop'])

        if random.random() < 0.1:
            if ga_par.mutation['type'] == 'random':
                self.ga_mutation_random()
        else:
            #if sum(self.ep.data['error_output'][0]) > 0:
            self.ga_mutiation_SGO_error()
            #else:
            self.ga_mutiation_SGO_pressure()
            self.ga_mutiation_SGO_flow()
            self.ga_mutiation_SGO_end_tank_level()
            self.ga_mutiation_SGO_acceleration_pump_pressure()
            self.ga_mutiation_SGO_energy()

    def ga_mutation_random(self):

        low_int_val = ga_par.pop['int_range_low']
        high_int_val = ga_par.pop['int_range_high']
        low_float_val = ga_par.pop['float_range_low']
        high_float_val = ga_par.pop['float_range_high']

        specimen = copy.deepcopy(self.pop['pop'])

        for i in range(len(specimen)):
            if random.random() < ga_par.mutation['percent_probability'] / 100:
                species = copy.deepcopy(specimen[i])
                if ga_par.number['int_genes'] > 0:

                    num_genes = int(ga_par.number['int_genes']*\
                                    ga_par.mutation['percent_genes']/100)

                    genes_position =[]
                    for j in range(num_genes):
                        genes_position.append(random.randint(0, ga_par.number['int_genes']-1))

                    for ge_pos in genes_position:
                        species[1][ge_pos] = random.randint(low_int_val, high_int_val+1)

                    self.pop['mate_mut'].append(species)

                if ga_par.number['float_genes'] > 0:
                    if random.random() < ga_par.mutation['percent_probability']/100:

                       num_genes = int(
                           ga_par.number['float_genes'] * ga_par.mutation['percent_genes'] / 100)

                       genes_position = []
                       for j in range(num_genes):
                           genes_position.append(random.randint(0, ga_par.number['float_genes']-1))

                       for ge_pos in genes_position:
                           if random.random()<0.5:
                                species[2][ge_pos] += random.random() * ga_par.mutation['random_delta']
                           else:
                                species[2][ge_pos] -= random.random() * ga_par.mutation['random_delta']

                           if species[2][ge_pos] > high_float_val:
                               species[2][ge_pos] = high_float_val
                           if species[2][ge_pos] < low_float_val:
                               species[2][ge_pos] = low_float_val
                    self.pop['mate_mut'].append(self.set_specimen_bound(species))

    def ga_mutiation_SGO_error(self):

        time_dur = sw_par.time['duration_h']
        for i in range(len(self.specimen)):
            if random.random() < ga_par.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if ga_par.number['float_genes'] > 0:# and species[8] == 0:
                    for t in range(sw_par.number['hydraulic_steps']):
                        error_value = species[3][t]
                        if error_value == 4:
                            for k in range(sw_par.number['pumps']):
                                flow_pump = species[6][k][t]
                                if flow_pump < 0.1:
                                    species[2][t + k * time_dur] += random.random() * ga_par.specialize_operators['error_delta_value']

                        species[8] = 1

                self.pop['mate_mut'].append(self.set_specimen_bound(species))

    def ga_mutiation_SGO_pressure(self):

        time_dur = sw_par.time['duration_h']
        for i in range(len(self.specimen)):
            if random.random() < ga_par.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if ga_par.number['float_genes'] > 0:# and species[8] == 0:
                    for j in range(time_dur):
                        kk=0
                        for k in range(len(sw_par.mes['nodes_names'])):

                            pressure_min_lev = sw_par.level['min_pressure'][k]
                            pressure_max_level = sw_par.level['max_pressure'][k]

                            if k < len(sw_par.mes['nodes_names']) - len(sw_par.tanks['names']):
                                pressure_WMO = species[4][k][j]
                                delta_min_pressure = - pressure_min_lev + pressure_WMO
                                delta_max_pressure = - pressure_WMO + pressure_max_level
                            else:
                                tank_level = species[5][kk][j]
                                delta_min_pressure = - pressure_min_lev + tank_level
                                delta_max_pressure = - tank_level + pressure_max_level
                                kk = kk + 1

                            if k == 0:
                                if delta_min_pressure < 0:
                                    species[2][j] += random.random() * ga_par.specialize_operators[
                                    'pressure_delta_value']
                                if delta_max_pressure < 0:
                                    species[2][j] -= random.random() * ga_par.specialize_operators[
                                        'pressure_delta_value']

                            if k == 1:
                                if  delta_min_pressure < 0:
                                    species[2][j] += random.random() * ga_par.specialize_operators[
                                    'pressure_delta_value']

                                if delta_max_pressure < 0:
                                    species[2][j] -= random.random() * ga_par.specialize_operators[
                                        'pressure_delta_value']

                             #Zbiorniki
                            if k == 2:
                                if delta_min_pressure < 0:

                                    for jj in range(1, ga_par.specialize_operators['bound_tank_level_horizon']+1):
                                        species[2][j - jj] += random.random() * \
                                                         ga_par.specialize_operators['pressure_delta_value']

                                if delta_max_pressure < 0:
                                    for jj in range(1, ga_par.specialize_operators['bound_tank_level_horizon']+1):
                                        species[2][j - jj] -= random.random() * \
                                                         ga_par.specialize_operators['pressure_delta_value']

                            species[8] = 1

                self.pop['mate_mut'].append(self.set_specimen_bound(species))

    def ga_mutiation_SGO_acceleration_pump_pressure(self):

        time_dur = sw_par.time['duration_h']
        for i in range(len(self.specimen)):
            if random.random() < ga_par.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if ga_par.number['float_genes'] > 0:# and species[8] == 0:
                    for j in range(time_dur-1):
                        kk = 0
                        for k in sw_par.mes_pressure_index['pump']:

                            max_delta_pump_pressure = ga_par.specialize_operators['max_delta_pump_pressure']

                            delta_pump_pressure = species[4][k][j] - species[4][k][j+1]

                            if k == 0:

                                if delta_pump_pressure < - max_delta_pump_pressure:
                                    species[2][j] += random.random() * ga_par.specialize_operators[
                                                                        'pressure_delta_value']
                                    species[2][j+1] -= random.random() * ga_par.specialize_operators[
                                                                        'pressure_delta_value']
                                if delta_pump_pressure > max_delta_pump_pressure:
                                    species[2][j] -= random.random() * ga_par.specialize_operators[
                                                                        'pressure_delta_value']
                                    species[2][j] += random.random() * ga_par.specialize_operators[
                                                                        'pressure_delta_value']

                            species[8] = 1

                self.pop['mate_mut'].append(self.set_specimen_bound(species))

    def ga_mutiation_SGO_end_tank_level(self):

        time_dur = sw_par.time['duration_h']
        for i in range(len(self.specimen)):
            if random.random() < ga_par.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if ga_par.number['float_genes'] > 0:#  and species[8] == 0:

                    for j in range(len(sw_par.tanks['names'])):

                        init_level = species[5][j][0]
                        end_level = species[5][j][-1]
                        if end_level - init_level > ga_par.specialize_operators['delta_initial_end_tank_level']:
                            for k in range(time_dur - ga_par.specialize_operators['end_tank_level_horizon'],
                                           time_dur):
                                species[2][j + k] -= random.random() * ga_par.specialize_operators[
                                    'pressure_delta_value']


                        if end_level - init_level < - ga_par.specialize_operators['delta_initial_end_tank_level']:
                            for k in range(time_dur - ga_par.specialize_operators['end_tank_level_horizon'],
                                           time_dur):
                                species[2][j + k] += random.random() * ga_par.specialize_operators[
                                    'pressure_delta_value']

                        species[8] = 1

                self.pop['mate_mut'].append(self.set_specimen_bound(species))

    def ga_mutiation_SGO_flow(self):

        time_dur = sw_par.time['duration_h']
        for i in range(len(self.specimen)):
            if random.random() < ga_par.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if ga_par.number['float_genes'] > 0:# and species[8] == 0:
                    for j in range(time_dur):
                        for k in range(len(sw_par.mes['links_names'])):

                            flow_mes = species[4][k][j]

                            flow_min_level = sw_par.level['min_flow'][k]
                            flow_max_level = sw_par.level['max_flow'][k]

                            delta_min_flow = - flow_min_level + flow_mes
                            delta_max_flow = - flow_mes + flow_max_level
                            if k == 0:
                                if delta_min_flow < 0:
                                    species[2][j] += random.random() * ga_par.specialize_operators[
                                    'flow_delta_value']

                                if delta_max_flow < 0:
                                    species[2][j] -= random.random() * ga_par.specialize_operators[
                                        'flow_delta_value']

                            species[8] = 1

                self.pop['mate_mut'].append(self.set_specimen_bound(species))

    def ga_mutiation_SGO_energy(self):

        time_dur = sw_par.time['duration_h']
        for i in range(len(self.specimen)):
            if random.random() < ga_par.mutation['percent_probability'] / 100:
                species = copy.deepcopy(self.specimen[i])
                if ga_par.number['float_genes'] > 0:# and species[8] == 0:
                    for j in range(time_dur):
                            if sw_par.energy_taryf[j] == sw_par.lc:
                                species[2][j] += random.random() * ga_par.specialize_operators[
                                    'energy_delta_value']

                            if sw_par.energy_taryf[j] == sw_par.hc:
                                species[2][j] -= random.random() * ga_par.specialize_operators[
                                    'energy_delta_value']

                            species[8] = 1

                self.pop['mate_mut'].append(self.set_specimen_bound(species))

    def set_specimen_bound(self, species):
        low_float_val = ga_par.pop['float_range_low']
        high_float_val = ga_par.pop['float_range_high']

        for t in range(len(species[2])):
            if species[2][t] > high_float_val:
                species[2][t] = high_float_val

            if species[2][t] < low_float_val:
                species[2][t] = low_float_val

        return species

    def ga_crossover(self):


        temp_cross = copy.deepcopy(self.pop['pop'])

        for te in self.pop['mate_mut']:
            temp_cross.append(te)

        self.pop['mate_cros'] = []
        if ga_par.crossover['type'] == 'single_point':
            for i in range(ga_par.number['specimen']):
                if random.random() < ga_par.crossover['percent_probability'] / 100:

                    specimen_to_crrossover = [random.randint(1, len(temp_cross)-1)]
                    specimen_to_crrossover.append(random.randint(1, len(temp_cross)-1))

                    specimen_0 = list(copy.deepcopy(temp_cross[specimen_to_crrossover[0]]))
                    specimen_1 = list(copy.deepcopy(temp_cross[specimen_to_crrossover[1]]))

                    if ga_par.number['int_genes'] > 0:
                        cros_position = random.randint(1, ga_par.number['int_genes'])

                        first_species_first_half = specimen_0[1][:cros_position[0]]
                        first_species_second_half = specimen_0[1][cros_position[0]:]
                        second_species_first_half = specimen_1[1][:cros_position[0]]
                        second_species_second_half = specimen_1[1][cros_position[0]:]

                        specimen_0[1] = second_species_first_half + first_species_second_half
                        specimen_1[1] = first_species_first_half + second_species_second_half

                    if ga_par.number['float_genes'] > 0:
                        cros_position = [random.randint(1, ga_par.number['float_genes'])]

                        first_species_first_half = specimen_0[2][:cros_position[0]]
                        first_species_second_half = specimen_0[2][cros_position[0]:]
                        second_species_first_half = specimen_1[2][:cros_position[0]]
                        second_species_second_half = specimen_1[2][cros_position[0]:]

                        specimen_0[2] = second_species_first_half + first_species_second_half
                        specimen_1[2] = first_species_first_half + second_species_second_half

                    self.pop['mate_cros'].append(specimen_0)
                    self.pop['mate_cros'].append(specimen_1)


        if ga_par.crossover['type'] == 'two_points':
            for i in range(ga_par.number['specimen']):
                if random.random() < ga_par.crossover['percent_probability'] / 100:

                    specimen_to_mutation = [random.randint(0, len(temp_cross)-1)]
                    specimen_to_mutation.append(random.randint(0, len(temp_cross)-1))

                    specimen_0 = list(copy.deepcopy(temp_cross[specimen_to_mutation[0]]))
                    specimen_1 = list(copy.deepcopy(temp_cross[specimen_to_mutation[1]]))

                    if ga_par.number['int_genes'] > 0:
                        cros_position = random.randint(1, ga_par.number['int_genes'])

                        first_species_first_part = specimen_0[1][:cros_position[0]]
                        first_species_middle_part = specimen_0[1][cros_position[0]:cros_position[1]]
                        first_species_last_part = specimen_0[1][cros_position[1]:]
                        second_species_first_part = specimen_0[1][:cros_position[0]]
                        second_species_middle_part = specimen_0[1][cros_position[0]:cros_position[1]]
                        second_species_last_part = specimen_0[1][cros_position[1]:]

                        specimen_0[1] = first_species_first_part + second_species_middle_part + first_species_last_part
                        specimen_1[1] = second_species_first_part + first_species_middle_part + second_species_last_part

                    if ga_par.number['float_genes'] > 0:
                        cros_position_1 = random.randint(1, ga_par.number['float_genes'])
                        cros_position_2 = random.randint(cros_position_1, ga_par.number['float_genes'])

                        first_species_first_part = specimen_0[1][:cros_position_1]
                        first_species_middle_part = specimen_0[1][cros_position_1:cros_position_2]
                        first_species_last_part = specimen_0[1][cros_position_2:]
                        second_species_first_part = specimen_0[1][:cros_position_1]
                        second_species_middle_part = specimen_0[1][cros_position_1:cros_position_2]
                        second_species_last_part = specimen_0[1][cros_position_2:]

                        specimen_0[1] = first_species_first_part + second_species_middle_part + first_species_last_part
                        specimen_1[1] = second_species_first_part + first_species_middle_part + second_species_last_part

                    self.pop['mate_cros'].append(specimen_0)
                    self.pop['mate_cros'].append(specimen_1)

        return 1

    def add_data_to_pop(self,epa,population_name,specimen_number):

        #3
        self.pop[population_name][specimen_number][3] = epa.data['error_output'][0]

        temp_pressure = []#4
        for i in range(len(epa.sw_parameters.mes['nodes_names']) - len(epa.sw_parameters.tanks['names'])):
            temp_pressure.append(epa.data['pressure_output_' + epa.sw_parameters.mes['nodes_names'][i]][0])
        self.pop[population_name][specimen_number][4] = temp_pressure

        temp_tank = []#5
        for i in range(len(epa.sw_parameters.tanks['names'])):
            temp_tank.append(epa.data['tank_output_' + epa.sw_parameters.tanks['names'][i]][0])

        self.pop[population_name][specimen_number][5] = temp_tank

        temp_flow = []#6
        for ln in epa.sw_parameters.mes['links_names']:
            temp_flow.append(epa.data['flow_output_' + ln][0])
        self.pop[population_name][specimen_number][6] = temp_flow

        #6 - Not repeted SGO mutation
        self.pop[population_name][specimen_number][8] = 0

    def ga_fitnes_function(self):

        if self.pop['num_iteration'] == 0:
            keys_list = ['pop']
        else:
            keys_list = ['pop', 'mate_mut', 'mate_cros']

        pomoc =[]
        for k in range(len(keys_list)):
            for i in (range(len(self.pop[keys_list[k]]))):
                species = self.pop[keys_list[k]][i][2]

                self.ep.get_data(species, ga_par.move_time)

                self.add_data_to_pop(self.ep, keys_list[k], i)

                penalties = {}
                penalties['error'] = ga_par.penalty_function_weights[0][0] * \
                                     self.penalty_function_error(self.ep)
                penalties['energy'] = ga_par.penalty_function_weights[0][1] * \
                                      self.penalty_function_energy(self.ep)
                penalties['pressure'] = ga_par.penalty_function_weights[0][2] * \
                                    self.penalty_function_pressure_limits(self.ep)
                penalties['flow'] = ga_par.penalty_function_weights[0][3] * \
                                    self.penalty_function_flow_limits(self.ep)
                penalties['tank'] = ga_par.penalty_function_weights[0][4] * \
                                    self.penalty_function_tank_limits(self.ep)
                penalties['tank_final'] = ga_par.penalty_function_weights[0][5] * \
                                          self.penalty_function_tank_final_state(self.ep)
                penalties['pump_pressure_acc'] = ga_par.penalty_function_weights[0][6] * \
                                        self.penalty_function_acceletarion_pump_pressure(self.ep)
                penalties['pump_flow_acc'] = ga_par.penalty_function_weights[0][7] * \
                                        self.penalty_function_acceletarion_pump_flow(self.ep)
                penalties['pump_speed_acc'] = ga_par.penalty_function_weights[0][8] * \
                                        self.penalty_function_acceletarion_pump_speed(self.ep)

                fitnes_fun =  penalties['error'] +  penalties['energy'] + penalties['pressure'] + penalties['flow'] + \
                              penalties['tank'] + penalties['tank_final'] + penalties['pump_pressure_acc']+\
                              penalties['pump_flow_acc'] + penalties['pump_speed_acc']

                self.pop[keys_list[k]][i][0] = fitnes_fun
                self.pop[keys_list[k]][i][7] = penalties

        return self.pop

    def ga_elitism(self):
        if ga_par.elitism[0]['type'] == 'first_n':
            keys_list = ['pop', 'mate_mut', 'mate_cros']
            ii=0
            pop_list=[]
            pomoc=[]
            for key in keys_list:
                for i in range(len(self.pop[key])):
                    ii+=1
                    pomoc.append([self.pop[key][i][0], key])
                    pop_list.append(self.pop[key][i])
            self.pop['pop'] = sorted(pop_list)[:int(ga_par.number['specimen']*\
                                                    ga_par.elitism[0]['n_type_precent'])]

    def ga_selection(self):

        temp_selection = copy.deepcopy(self.pop['pop'])

        self.pop['selection'] = []

        if ga_par.selection[0]['type'] == 'tournament':
            for i in range(ga_par.number['specimen']):
                if random.random() < ga_par.selection[0]['precent_probability']:

                    specimen_to_tournament = [random.randint(0, ga_par.number['specimen']-1)]
                    specimen_to_tournament.append(random.randint(0, ga_par.number['specimen']-1))

                    specimen_0 = copy.deepcopy(temp_selection[specimen_to_tournament[0]])
                    specimen_1 = copy.deepcopy(temp_selection[specimen_to_tournament[1]])

                    if specimen_0[0] <= specimen_1[0]:
                        self.pop['selection'].append(specimen_0)

    def ga_statistics(self):

        self.pop['last_best_solution'] = self.pop['pop'][0]

        self.pop['best_solutions'].append(self.pop['last_best_solution'][0])

        population_fun_fitnes = []

        for i in range(ga_par.number['specimen']):
            population_fun_fitnes.append(self.pop['pop'][i][0])

        self.pop['mean_pop_value'].append(mean(population_fun_fitnes))

    def ga_stop_criterion(self):

        output = 1
        if self.pop['num_iteration'] == 0:
            self.last_solution = [self.pop['last_best_solution'][0]]
            self.diff_fitness_function = [0]

        else:
            self.last_solution.append(self.pop['last_best_solution'][0])
            self.diff_fitness_function.append(abs(self.last_solution[-1] - self.last_solution[-2]))

            last_n_diff = self.diff_fitness_function[-ga_par.stop_criterion['min_generations']:]

        if self.pop['num_iteration'] >= ga_par.stop_criterion['min_generations'] and\
            self.pop['num_iteration'] >0:

            if sum(last_n_diff) < ga_par.stop_criterion['min_change'] and \
                sum(self.pop['last_best_solution'][3]) == 0:

                output = 0

        self.pop['num_iteration'] += 1

        return output

    def print_statistics(self):

        self.ga_statistics()
        self.print_stats = f"It. num.: {ga.pop['num_iteration']}, "+\
                          f"Mean: {ga.pop['mean_pop_value'][-1]:.0f}, "+ \
                          f"Best: {ga.pop['last_best_solution'][0]:.0f}, "+\
                          f"Energy: {ga.pop['last_best_solution'][7]['energy']:.2f}, "+\
                          f"Pres: {ga.pop['last_best_solution'][7]['pressure']:.0f}, "+\
                          f"Flow: {ga.pop['last_best_solution'][7]['flow']:.0f}, "+\
                          f"Tank: {ga.pop['last_best_solution'][7]['tank']:.0f}, "+\
                          f"Tank final: {ga.pop['last_best_solution'][7]['tank_final']:.0f}, " + \
                          f"Pres\Flow\Speed acc: {ga.pop['last_best_solution'][7]['pump_pressure_acc']:.0f}, " + \
                          f"{ga.pop['last_best_solution'][7]['pump_flow_acc']:.0f}, " + \
                          f"{ga.pop['last_best_solution'][7]['pump_speed_acc']:.0f}"
        print('--------------------------------------------------------------------------------------------')

        print(self.print_stats)

    def plot_mes(self, axis, type, which, location):

        x = np.arange(0, sw_par.time['duration_h']*int(3600/sw_par.time['hydraulic_step_s']))

        if type == 4:
            mes_name = 'nodes'
            name = 'pressure'
            ylabel = 'Pressure [m]'
        if type ==5:
            mes_name = 'tank'
            name = 'tank'
            ylabel = 'Level [m]'
        if type == 6:
            mes_name = 'links'
            name = 'flow'
            ylabel = 'Flow [m^3/h]'

        upper_bound = [sw_par.level['max_'+name][which[0]]]*sw_par.time['duration_h']\
                      *int(3600/sw_par.time['hydraulic_step_s'])
        lower_bound = [sw_par.level['min_'+name][which[0]]]*sw_par.time['duration_h']\
                      *int(3600/sw_par.time['hydraulic_step_s'])
        self.ax[axis[0], axis[1]].clear()
        color = ['b','m','g','k']
        for i in range(len(which)):
            self.ax[axis[0], axis[1]].plot(x, self.pop['last_best_solution'][type][which[i]],
                                            color=color[i], linewidth=1)
        self.ax[axis[0], axis[1]].plot(x, lower_bound, color='r', linewidth=1)
        self.ax[axis[0], axis[1]].plot(x, upper_bound, color='r', linewidth=1)
        self.ax[axis[0], axis[1]].set_ylabel(ylabel, fontsize=8)
        self.ax[axis[0], axis[1]].set_xlabel("Time [h]", fontsize=8)


        labels =[]
        names =''
        for i in range(len(which)):
            labels.append(sw_par.mes[mes_name + '_names'][which[i]])
            if i < len(which)-1:
                names += sw_par.mes[mes_name +'_names'][which[i]]+', '
            else:
                names += sw_par.mes[mes_name +'_names'][which[i]]

        #self.ax[axis[0], axis[1]].set_title(name + ' at '+ names,
        #                                    fontsize = 8)
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
        figure, ax = plt.subplots(nrows=3, ncols=6, figsize = (18,10))
        mngr = plt.get_current_fig_manager()
        #mngr.window.geometry("+100+0")
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
            if er > 0:
                error += str(i) +'/'+ str(er)+', '
            i += 1
        self.figure.text(0.01, 0.95, ' ' * 800, fontsize=12, bbox=dict(facecolor='white', edgecolor='white'))
        self.figure.text(0.01, 0.95, self.print_stats, fontsize=12, bbox ={'facecolor':'white'})
        self.figure.text(0.01, 0.92, ' '*100, fontsize=12, bbox=dict(facecolor= 'white', edgecolor = 'white'))
        self.figure.text(0.01, 0.92, error, fontsize=12, bbox={'facecolor': 'white'})

    def ga_plot(self,i, epa):

        self.plot_fitnes_values([0,0],i)
        self.plot_mes([0, 1], 4, [0],'upper right')
        self.plot_mes([0, 2], 4, [1], 'upper right')
        self.plot_mes([0, 3], 4, [2], 'upper right')
        self.plot_mes([0, 4], 4, [3], 'upper right')
        self.plot_mes([0, 5], 4, [4], 'upper right')
        self.plot_mes([1, 0], 4, [5], 'upper right')
        self.plot_mes([1, 1], 4, [6], 'upper right')
        self.plot_mes([1, 2], 4, [7], 'upper right')
        self.plot_mes([1, 3], 4, [8], 'upper right')
        self.plot_mes([1, 4], 4, [9], 'upper right')
        self.plot_mes([1, 5], 5, [0], 'upper right')
        self.plot_mes([2, 0], 5, [1,2], 'upper right')
        self.plot_mes([2, 2], 5, [3], 'upper right')
        self.plot_mes([2, 3], 6, [0, 1, 3], 'upper right')
        self.plot_mes([2, 4], 6, [4,5], 'upper right')
        self.plot_mes([2, 5], 6, [5], 'upper right')
        #self.plot_mes([1, 0], 2, [0], 'upper right')
        self.print_values_on_plot(epa)

        self.save_figure()

        self.figure.canvas.draw()

        self.figure.canvas.flush_events()

        time.sleep(0.1)

    def save_figure(self):
        self.figure.savefig('last_figure.png')

def save_variabele_space_to_file(file_name, input_dictionary):

    key = 'pop'
    temp_date = []
    for i in range(ga_par.number['specimen']):
        temp_date.append(f"{input_dictionary['pop'][i][2:]}\n")

    with open(file_name,'a') as f:
        f.writelines(temp_date)
        f.close()

if __name__ == '__main__':
        print('\nTo jest biblioteka pomocna przy algorytmie genetycznym')

        ga = GA_function()

        ga.ga_initialization()


        ga.ga_fitnes_function()

        ga.ga_statistics()

        ga.ga_init_plot()

        for i in range(1,ga_par.number['generations']+1):

            #ga.ga_selection()

            ga.ga_mutation()

            ga.ga_crossover()

            ga.ga_fitnes_function()

            ga.ga_elitism()

            #save_variabele_space_to_file('last_best_pop_'+ga_par.data + '.txt', ga.pop)

            ga.print_statistics()

            ga.ga_plot(i, ga.pop)

            ga.ep.save_last_net(ga.pop['pop'][0][2])

            if ga.ga_stop_criterion() == 0:

                break