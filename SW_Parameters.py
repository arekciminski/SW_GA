from Get_data_from_DWDS import Epa

class sw_par:

    file_name = "enag_model.inp",

    tanks = {'names': ['4', '5'],
             'initial_lev': [6, 5]}

    lc = 0.84
    hc = 2*lc
    energy_taryf = [lc]*6
    energy_taryf.extend([hc]*4)
    energy_taryf.extend([lc] * 4)
    energy_taryf.extend([hc] * 4)
    energy_taryf.extend([lc] * 6)

    mes = {'links_names': ['110', '210','323'],
            'nodes_names': ['10','23','32', '35', '36', '4','5']}

    mes['tanks_names'] = tanks['names']

    mes_pressure_index = {'pump': [0,1],
                      'tank': [5,6],
                      'inside': [2,3,4]}

    pumps = {'names': ['110','210','323'],
             'patterns_names': ['2','3','4']}

    demand = {'patterns_names': ['1']}

    time = {'duration_h': 24,  # [h]
            'hydraulic_step_s': 3600}

    time['duration_s'] = time['duration_h'] * time['hydraulic_step_s']  # [s]

    hydraulic_monitor_values = ['flow', 'energy', 'pressure']

    level = {'min_pressure': [40, 50, 40, 40, 40, 2, 2],
             'max_pressure': [55, 60, 60, 60, 60, 9, 19],
             'min_tank' : [2, 2],
             'max_tank': [9, 19],
             'min_flow': [0,0,0],
             'max_flow': [200, 200, 200],
             'max_pumps_speed':[1],
             'min_pumps_speed': [0]
    }

    number = {'mes_links' : len(mes['links_names']),
            'mes_nodes' : len(mes['nodes_names']),
            'pumps' : len(pumps['names']),
            'demands' : len(demand['patterns_names']),
            'tanks' : len(tanks['names'])}

