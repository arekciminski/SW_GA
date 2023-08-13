
class sw_par:

    file_name = "chojnice_kwiecien_obiekt.inp",

    tanks = {'names': ['180'],
             'initial_lev': [2.58]}

    lc = 0.84
    hc = 2*lc
    energy_taryf = [lc]*6
    energy_taryf.extend([hc]*4)
    energy_taryf.extend([lc] * 4)
    energy_taryf.extend([hc] * 4)
    energy_taryf.extend([lc] * 6)

    mes = {'links_names': ['F1', 'K1', 'P1'],
            'nodes_names': ['001', '002', '055', '083', '158', '180']}
    mes['tank_names'] = tanks['names']

    mes_head_index = {'pump': [0, 1, 2],
                      'tank': [5],
                      'inside': [3, 4]}

    pumps = {'names': ['F1', 'K1', 'P1'],
             'patterns_names': ['f1', 'k1', 'p1']}

    demand = {'patterns_names': ['par2', 'kar2']}

    time = {'duration_h': 24,  # [h]
            'hydraulic_step_s': 3600}

    time['duration_s'] = time['duration_h'] * time['hydraulic_step_s']  # [s]

    hydraulic_monitor_values = ['flow', 'energy', 'head']

    level = {'min_head': [150, 150, 150, 196, 195, 167.2],
             'max_head': [185, 205, 205, 205, 205, 170.2],
             'min_tank' : [167.2],
             'max_tank': [170.8],
             'min_flow': [0, 0, 0],
             'max_flow': [250, 250, 140]
    }

    number = {'mes_links' : len(mes['links_names']),
            'mes_nodes' : len(mes['nodes_names']),
            'pumps' : len(pumps['names']),
            'demands' : len(demand['patterns_names']),
            'tanks' : len(tanks['names'])}
