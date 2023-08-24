class sw_par:

    file_name = "Gdynia_Dolna_strefa_probny",

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

    tanks = {'names': ['ZB_Witomino', 'ZB_Kwidzynski_L','ZB_Kwidzynski_P','Zb_P1'],
             'initial_lev': [2, 2.07, 2.07, 2.5]}

    pumps = {'names': ['P1_Sieradzka','SUW_Kolibki', 'P_P1','SUW_REDA']}
    pumps['patterns_names'] = pumps['names']

    mes = {'links_names': pumps['names'],
            'nodes_names': ['M33_3', '309','443','M38_1','316', '312', '373', 'R1308_1',
                            '303-pob', 'M35_Partyz_1pob', 'M89_1-pob','372-pob','M76_16','RUMIAREDA_POBOR','1']}
    mes['nodes_names'].extend(tanks['names'])
    mes['tank_names'] = tanks['names']

    level = {'min_pressure': [40, 46, 42, 42, 20, 30, 28, 20,
                              64, 58, 55, 64, 56, 45, 40],
             'max_pressure': [52, 54, 49, 56, 35, 40, 34, 32,
                              72, 66, 80, 70, 64, 55, 55],
             'min_tank' : [1, 1, 1, 1],
             'max_tank': [4.0, 4.25, 4.25, 5.0],
             'min_flow': [0, 0, 0, 0],
             'max_flow': [4000, 3000, 17000, 30000],
             'max_pumps_speed':[1.4],
             'min_pumps_speed': [0]
    }

    level['min_pressure'].extend(level['min_tank'])
    level['max_pressure'].extend(level['max_tank'])

    mes_pressure_number = {'pump': range(0,len(pumps['names'])),
                          'inside': range(len(pumps['names']),len(mes['nodes_names']) - len(tanks['names'])),
                          'tank': range(len(mes['nodes_names']) - len(tanks['names']),len(mes['nodes_names']))}



    hydraulic_monitor_values = ['flow', 'energy', 'pressure']

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