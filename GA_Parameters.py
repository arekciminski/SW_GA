from SW_Parameters import sw_par
import datetime
class ga_par:
    data = str(datetime.datetime.now()).split('.')
    data = data[0].replace(':', '_').replace(' ', '_')

    #[error, energy, pressure, flow, tank, tank_initial,
    # accelaration pump pressure, accelaration pump flow,
    # accelaration pump speed]
    penalty_function_weights = [1e14,
                                10,
                                50,
                                10,
                                10,
                                10,
                                0,
                                0,
                                0],
    number = {'generations': 500,
              'specimen': 10,
              'float_genes': int(len(sw_par.pumps['names']) * sw_par.time['duration_h'] * \
                             3600/sw_par.time['hydraulic_step_s']),
              'int_genes': 0}


    elitism ={ 'type': 'first_n',
               'n_type_precent': 1},

    selection = {'type': 'tournament',
                 'precent_probability': 0.5},  # 'first_n',roulete  # rws, sus, random, tournament

    crossover = {'percent_probability': 20,
                 'type': 'single_point'}  #single_point, two_points, uniform , scattered
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
    mutation = {'percent_genes': 10,
                'percent_probability': 90,
                'type': 'random',  # random, inversion , scramble ,
                'random_delta': 0.15}

    pop = {'float_range_low': 0,
           'float_range_high': 1.3,
           'int_range_low': 0,
           'int_range_high': 1,
           }

    stop_criterion = {
        'min_generations': 20,
        'min_change': 1e-1}

    move_time = 0