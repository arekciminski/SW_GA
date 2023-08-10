from SW_Parameters import sw_par

class ga_par:
    #[error, energy, head, flow, tank, tank_initial]
    penalty_function_weights = [1e9,
                                100,
                                1000,
                                100,
                                10000       ,
                                1000],
    number = {'generations': 500,
           'specimen': 100,
           'float_genes': len(sw_par.pumps['names']) * sw_par.time['duration_h'],
            'int_genes': 0}


    selection ={ 'type': 'first_n'},#'first_n',roulete  # rws, sus, random, tournament
    crossover = {'percent_probability': 30,
                 'type': 'single_point'}  # two_points , uniform , scattered
    specialize_operators ={'error_delta_value': 0.25,
                            'head_delta_value': 0.2,
                            'flow_delta_value': 0.2}
    mutation = {'percent_genes': 20,
                'percent_probability': 90,
                'type': 'random',  # random, inversion , scramble ,
                'random_delta': 0.05}

    pop = {'float_range_low': 0.5,
            'float_range_high': 1,
           'int_range_low': 0,
           'int_range_high': 1,
           }

    stop_criterion = {
        'min_generations': 10,
        'min_change': 1e-2}