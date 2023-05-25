from epanettools import epanet2 as et

class Epa:

    def __init__(self, SW_parameters):
        self.SW_parameters = SW_parameters

    def open_epanet(self):
        ret = et.ENopen(self.SW_parameters['file_name'], "Net3.rpt", "")

    def close_epanet(self):
        ret = et.ENclose()

    def get_number_of_nodes(self):
        ret, nnodes = et.ENgetcount(et.EN_NODECOUNT)
        self.nnodes = nnodes
        return nnodes

    def get_number_of_links(self):
        ret, nlinks = et.ENgetcount(et.EN_LINKCOUNT)
        self.nlinks = nlinks
        return nlinks

    def get_link_index(self):
        link_index = []
        for ln in self.SW_parameters['mes_links_names']:
            ret, lx = et.ENgetlinkindex(ln)
            link_index.append(lx)
        self.link_index = link_index

        pump_index = []
        for ln in self.SW_parameters['pumps_names']:
            ret, lx = et.ENgetlinkindex(ln)
            pump_index.append(lx)
        self.pumps_index = pump_index

    def get_node_index(self):
        node_index = []
        for nn in self.SW_parameters['mes_nodes_names']:
            ret, nx = et.ENgetnodeindex(nn)
            node_index.append(nx)
            self.node_index = node_index

        tank_index = []
        for nn in self.SW_parameters['tanks_names']:
            ret, nx = et.ENgetnodeindex(nn)
            tank_index.append(nx)
            self.tank_index = tank_index

    def get_pattern_index(self):
        pattern_index = []
        for pn in self.SW_parameters['pumps_patterns_names']:
            ret, px = et.ENgetpatternindex(pn)
            pattern_index.append(px)
        self.pumps_pattern_index = pattern_index

        pattern_index = []
        for pn in self.SW_parameters['demand_patterns_names']:
            ret, px = et.ENgetpatternindex(pn)
            pattern_index.append(px)
        self.demand_pattern_index = pattern_index

    def set_time_duration(self):
        et.ENsettimeparam(0, self.SW_parameters['time_duration_s'])

    def set_tank_inital(self):
        for tank_in, initial_val in zip(self.tank_index, self.SW_parameters['initial_tanks_lev']):
            et.ENsetnodevalue(tank_in, 8, initial_val)

    def set_patern_values(self,species):
        for i in range(self.SW_parameters['num_pumps']):
             for j in range(self.SW_parameters['time_duration_h']):
                 et.ENsetpatternvalue(self.pumps_pattern_index[i], j + 1, species[i*\
                                                                        self.SW_parameters['time_duration_h'] + j])

    ''' for i in range(self.parameters['num_demands']):
            for j in range(self.parameters['time_duration_h']):
                et.ENsetpatternvalue(self.demand_pattern_index[i], j + 1, self.random['demand_val'][i][j])'''

    def save_temp_file(self):
        et.ENsaveinpfile('temporary_' + self.SW_parameters['file_name'])

    def get_set_parameters(self):
        self.get_link_index()
        self.get_node_index()
        self.get_pattern_index()
        self.set_time_duration()
        #self.save_temp_file()

    def prepare_empty_dict_to_comput(self):

        self.data = {'error_output': []}

        self.data['time_output'] = []

        for mes_node in self.SW_parameters['mes_nodes_names']:
            self.data['head_output_' + mes_node] = []

        for mes_link in self.SW_parameters['mes_links_names']:
            self.data['flow_output_' + mes_link] = []

        for mes_tank in self.SW_parameters['tanks_names']:
            self.data['tank_output_' + mes_tank] = []

        for pumps_name in self.SW_parameters['pumps_names']:
            self.data['energy_output_' + pumps_name] = []

        for tank_input in self.SW_parameters['tanks_names']:
            self.data['tank_input_' + tank_input] = []

        for pump_input in self.SW_parameters['pumps_names']:
            self.data['pump_input_' + pump_input] = []

        for demand_input in self.SW_parameters['demand_patterns_names']:
            self.data['demand_input_' + demand_input] = []

    def get_hydraulic_values(self):

        et.ENopenH()
        et.ENinitH(0)
        time = []
        flow = [[] for i in range(self.SW_parameters['num_mes_links'])]
        head = [[] for i in range(len(self.node_index))]
        energy = [[] for i in range(self.SW_parameters['num_pumps'])]
        error = []
        while True:
            ret, t = et.ENrunH()
            if t % self.SW_parameters['hydraulic_step_s'] == 0:
                error.append(ret)
                time.append(t)
                if self.SW_parameters['hydraulic_monitor_values'][0] == 'flow':
                    for i in range(0, len(self.link_index)):
                        ret, p = et.ENgetlinkvalue(self.link_index[i], et.EN_FLOW)
                        flow[i].append(p)

                if self.SW_parameters['hydraulic_monitor_values'][2] == 'head':
                    for i in range(0, len(self.node_index)):
                        ret, p = et.ENgetnodevalue(self.node_index[i], et.EN_HEAD)
                        head[i].append(p)

                if self.SW_parameters['hydraulic_monitor_values'][1] == 'energy':
                    for i in range(0, self.SW_parameters['num_pumps']):
                        ret, p = et.ENgetlinkvalue(self.pumps_index[i], et.EN_ENERGY)
                        energy[i].append(p)

            ret, tstep = et.ENnextH()
            if (tstep <= 0):
                break
        ret = et.ENcloseH()
        return flow, energy, head, error, time

    def insert_data(self, head, flow, energy, error, time):

        for ii in range(0, self.SW_parameters['num_mes_nodes'] - self.SW_parameters['num_tanks'] + 1):
            self.data['head_output_' + self.SW_parameters['mes_nodes_names'][ii]].append(head[ii][:-1])

        for ii in range(self.SW_parameters['num_mes_nodes'] - self.SW_parameters['num_tanks'],
                        self.SW_parameters['num_mes_nodes']):
            self.data['tank_output_' + self.SW_parameters['tanks_names'][
                self.SW_parameters['num_mes_nodes'] - ii - 1]].append(head[ii][1:])

        for ii in range(0, self.SW_parameters['num_tanks']):
            temp_tank_init = []
            for i in range(0, self.SW_parameters['time_duration_h']):
                temp_tank_init.append(self.SW_parameters['initial_tanks_lev'][ii])
            self.data['tank_input_' + self.SW_parameters['tanks_names'][ii]].append(temp_tank_init)

        for ii in range(self.SW_parameters['num_mes_links']):
            self.data['flow_output_' + self.SW_parameters['mes_links_names'][ii]].append(flow[ii][:-1])

        for ii in range(self.SW_parameters['num_pumps']):
            self.data['energy_output_' + self.SW_parameters['pumps_names'][ii]].append(energy[ii][:-1])

        self.data['time_output'].append(time[:-1])

        self.data['error_output'].append(error[:-1])

    def get_data(self,species):

        self.prepare_empty_dict_to_comput()

        self.open_epanet()
        self.get_set_parameters()

        self.set_tank_inital()
        self.set_patern_values(species)
        [flow, energy, head, error, time] = self.get_hydraulic_values()

        self.insert_data(head, flow, energy, error, time)

        self.close_epanet()


if __name__ == '__main__':
    print('To jest biblioteka pomocna przy generowaniu wyników z Epanetu')
