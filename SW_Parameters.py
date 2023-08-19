from Get_data_from_DWDS import Epa

class sw_par:

    file_name = "siec_maly.inp",

    tanks = {'names': ['7'],
             'initial_lev': [2]}

    lc = 0.84
    hc = 2*lc
    energy_taryf = [lc]*6
    energy_taryf.extend([hc]*4)
    energy_taryf.extend([lc] * 4)
    energy_taryf.extend([hc] * 4)
    energy_taryf.extend([lc] * 6)

    mes = {'links_names': ['1'],
            'nodes_names': ['2', '6', '7']}
    mes['tank_names'] = tanks['names']

    mes_pressure_index = {'pump': [0],
                      'tank': [2],
                      'inside': [1]}

    pumps = {'names': ['1'],
             'patterns_names': ['pompka']}

    demand = {'patterns_names': ['demand']}

    time = {'duration_h': 24,  # [h]
            'hydraulic_step_s': 3600}

    time['duration_s'] = time['duration_h'] * time['hydraulic_step_s']  # [s]

    hydraulic_monitor_values = ['flow', 'energy', 'pressure']

    level = {'min_pressure': [30, 31, 0.5],
             'max_pressure': [35, 36, 4.5],
             'min_tank' : [0.5],
             'max_tank': [4.5],
             'min_flow': [0],
             'max_flow': [10],
             'max_pumps_speed':[1],
             'min_pumps_speed': [0]
    }

    number = {'mes_links' : len(mes['links_names']),
            'mes_nodes' : len(mes['nodes_names']),
            'pumps' : len(pumps['names']),
            'demands' : len(demand['patterns_names']),
            'tanks' : len(tanks['names'])}

