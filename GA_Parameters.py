from SW_Parameters import sw_par
import datetime
class ga_par:
    data = str(datetime.datetime.now()).split('.')
    data = data[0].replace(':', '_').replace(' ', '_')

    penalty_function_weights = [1e14, #error
                                1,    #energy
                                1000, #pressure
                                0,    #flow
                                10000,#tank
                                1000, #tank_inital
                                0.1, #accelaration pump pressure,
                                0,    #accelaration pump flow,
                                0],   #accelaration pump speed
    number = {'generations': 500,
              'specimen': 70,
              'float_genes': int(len(sw_par.pumps['names']) * sw_par.number['hydraulic_steps']),
              'int_genes': 0}


    elitism ={ 'type': 'first_n',
               'n_type_precent': 1},

    selection = {'type': 'tournament',
                 'precent_probability': 0.5},  # 'first_n',roulete  # rws, sus, random, tournament

    crossover = {'percent_probability': 20,
                 'type': 'intermediate',  #single_point, two_points, intermediate
                 'inremediate_d_par' : 0.5}

    specialize_operators ={'error_delta_value': 0.1,
                           'pressure_delta_value': 0.2,
                           'flow_delta_value': 0.2,
                           'end_tank_level_horizon': 8,
                           'bound_tank_level_horizon': 4,
                           'delta_initial_end_tank_level': 0.5,
                           'max_delta_pump_pressure': 3,
                           'max_delta_pump_flow': 5*sw_par.multiplication,
                           'max_delta_pump_speed': 0.3,
                           'energy_delta_value': 0.1}
    mutation = {'percent_genes': 10,
                'percent_probability': 90,
                'type': 'random',  # random, inversion , scramble ,
                'random_delta': 0.15}

    pop = {'float_range_low': 0,
           'float_range_high': 1.2,
           'int_range_low': 0,
           'int_range_high': 1,
           }

    show = {'number_last_solutions': 20}

    stop_criterion = {
        'min_generations': 5,
        'min_change': 1e-1}

    move_time = 0