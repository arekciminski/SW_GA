from SW_Parameters import sw_par

class ga_par:
    #[error, energy, head, flow, tank, tank_initial, accelaration pump head]
    penalty_function_weights = [1e9,
                                50,
                                1000,
                                100,
                                10000       ,
                                1000,
                                1],
    number = {'generations': 500,
           'specimen': 100,
           'float_genes': len(sw_par.pumps['names']) * sw_par.time['duration_h'],
            'int_genes': 0}


    selection ={ 'type': 'first_n',
                 'n_type_precent': 1},#'first_n',roulete  # rws, sus, random, tournament
    crossover = {'percent_probability': 40,
                 'type': 'single_point'}  # two_points , uniform , scattered
    specialize_operators ={'error_delta_value': 0.25,
                           'head_delta_value': 0.2,
                           'flow_delta_value': 0.2,
                           'end_tank_level_horizon': 5,
                           'delta_ini_end_tank_level': 0.5,
                           'max_delta_head': 10}
    mutation = {'percent_genes': 20,
                'percent_probability': 90,
                'type': 'random',  # random, inversion , scramble ,
                'random_delta': 0.15}

    pop = {'float_range_low': 0.5,
           'float_range_high': 1,
           'int_range_low': 0,
           'int_range_high': 1,
           }

    stop_criterion = {
        'min_generations': 10,
        'min_change': 1e-2}