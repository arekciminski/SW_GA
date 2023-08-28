class sw_par:

    file_name = "Gd_dec_pomp_hyd_chlor_read.inp",

    tanks = {'names': ['ZB_Witomino', 'ZB_Kwidzynski_L','ZB_Kwidzynski_P','Zb_P1'],
             'initial_lev': [2, 2.07, 2.07, 2.5]}

    time = {'duration_h': 24,  # [h]
            'hydraulic_step_s': 1800}
    multiplication = int(3600/time['hydraulic_step_s'])

    time['duration_s'] = time['duration_h'] * 3600  # [s]

    lc = 0.84
    hc = 2 * lc
    energy_taryf = [lc] * 6 * multiplication
    energy_taryf.extend([hc] * 4 * multiplication)
    energy_taryf.extend([lc] * 4 * multiplication)
    energy_taryf.extend([hc] * 4 * multiplication)
    energy_taryf.extend([lc] * 6 * multiplication)

    pumps = {'names': ['P1_SIERADZKA','SUW_KOLIBKI', 'P_P1','SUW_REDA']}
    pumps['patterns_names'] = pumps['names']

    mes = {'links_names': pumps['names'],
            'nodes_names': ['M33_3', '309', '443', 'M38_1', '316', '312', '373', 'R1308_1',
                            'M50_4', '37_2', 'M89_1-pob','372-pob','M28_1','RUMIAREDA_POBOR','1']}
    level = {'min_pressure': [30, 30, 30, 30, 20, 30, 28, 20, 1,  1,  53, 60, 10, 44, 40],
             'max_pressure': [52, 60, 49, 56, 35, 40, 34, 32, 6, 6, 80, 66, 25, 50, 55],
             'min_tank' : [1, 1, 1, 1],
             'max_tank': [5.0, 4.25, 4.25, 4.0],
             'min_flow': [0, 0, 0, 0, 0, 0, 0, 0],
             'max_flow': [3000, 3000, 20000, 32000, 600, 140, 200, 200],
             'max_pumps_speed':[1],
             'min_pumps_speed': [0]
    }


    mes['nodes_names'].extend(tanks['names'])
    mes['tank_names'] = tanks['names']

    mes_pressure_index = {'pump': range(0,len(pumps['names'])),
                          'inside': range(len(pumps['names']),len(mes['nodes_names']) - len(tanks['names'])),
                          'tank': range(len(mes['nodes_names']) - len(tanks['names']),len(mes['nodes_names']))}




    hydraulic_monitor_values = ['flow', 'energy', 'pressure']

    level['min_pressure'].extend(level['min_tank'])
    level['max_pressure'].extend(level['max_tank'])

    number = {'mes_links' : len(mes['links_names']),
              'mes_nodes' : len(mes['nodes_names']),
              'pumps' : len(pumps['names']),
              'tanks' : len(tanks['names']),
              'hydraulic_steps':int(time['duration_h'] * 3600 / time['hydraulic_step_s'])
              }

    if len(mes['nodes_names']) != len(level['min_pressure']):
        Warning('Number of names and bounds are not the same')

    if len(level['max_pressure']) != len(level['min_pressure']):
        Warning('Number of min and max bounds are not the same')