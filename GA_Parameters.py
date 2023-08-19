from SW_Parameters import sw_par
import datetime
class ga_par:
    data = str(datetime.datetime.now()).split('.')
    data = data[0].replace(':', '_').replace(' ', '_')

    #[error, energy, pressure, flow, tank, tank_initial,
    # accelaration pump pressure, accelaration pump flow,
    # accelaration pump speed]
    penalty_function_weights = [1e9,
                                10,
                                5000,
                                100,
                                100000,
                                10000,
                                10,
                                10,
                                0],
    number = {'generations': 500,
           'specimen': 50,
           'float_genes': len(sw_par.pumps['names']) * sw_par.time['duration_h'],
            'int_genes': 0}


    selection ={ 'type': 'first_n',
                 'n_type_precent': 1},#'first_n',roulete  # rws, sus, random, tournament
    crossover = {'percent_probability': 80,
                 'type': 'single_point'}  # two_points , uniform , scattered
    specialize_operators ={'error_delta_value': 0.3,
                           'pressure_delta_value': 0.3,
                           'flow_delta_value': 0.3,
                           'end_tank_level_horizon': 8,
                           'bound_tank_level_horizon': 4,
                           'delta_initial_end_tank_level': 0.5,
                           'max_delta_pump_pressure': 3,
                           'max_delta_pump_flow': 5,
                           'max_delta_pump_speed': 0.3,
                           'energy_delta_value': 0.2}
    mutation = {'percent_genes': 1,
                'percent_probability': 90,
                'type': 'random',  # random, inversion , scramble ,
                'random_delta': 0.15}

    pop = {'float_range_low': 0,
           'float_range_high': 1,
           'int_range_low': 0,
           'int_range_high': 1,
           }

    stop_criterion = {
        'min_generations': 10,
        'min_change': 1e-1}

    move_time = 0