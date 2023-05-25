import numpy as np

from GA_functions import GA_function

file_name = "chojnice_kwiecien_obiekt.inp"

SW_parameters = {
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


SW_parameters['num_mes_links'] = len(SW_parameters['mes_links_names'])
SW_parameters['num_mes_nodes'] = len(SW_parameters['mes_nodes_names'])
SW_parameters['num_pumps'] = len(SW_parameters['pumps_names'])
SW_parameters['num_demands'] = len(SW_parameters['demand_patterns_names'])
SW_parameters['num_tanks'] = len(SW_parameters['tanks_names'])
SW_parameters['time_duration_s'] = SW_parameters['time_duration_h']*3600 ##[s]

GA_parameters = {
    'penalty_function_weights': [1e5, 10, 20, 30, 40],
    'num_generations': 100,
    'num_parents_mating': SW_parameters['num_pumps'] * SW_parameters['time_duration_h'],
    'sol_per_pop': 50,
    'gene_type': float,
    'num_genes': SW_parameters['num_pumps'] * SW_parameters['time_duration_h'],
    'init_range_low': 0,
    'init_range_high': 1,
    'parent_selection_type': 'sss', #rws, sus, random, tournament
    'keep_parents': -1,
    'keep_elitism': 10,
    'crossover_probability': 0.9,
    'crossover_type': 'single_point', #two_points , uniform , scattered
    'mutation_type':  'random', #swap, inversion , scramble ,
    'mutation_percent_genes': 50,
    'random_mutation_min_val': 0,
    'random_mutation_max_val': 1,
    'on_mutation': None,
    'save_best_solutions': True,
    'parallel_processing': None,

}

species = {'pump_input': np.random.random_sample(size= (SW_parameters['num_pumps'], SW_parameters['time_duration_h']))}

ga = GA_function(GA_parameters,SW_parameters)

ga.fitnes_function(species)
