from epanettools import epanet2 as et
from SW_Parameters import sw_par

class Epa:

    def __init__(self, sw_parameters):
        self.sw_parameters = sw_parameters

    def open_epanet(self):
        ret = et.ENopen(self.sw_parameters.file_name[0], "Net3.rpt", "")

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
        for ln in self.sw_parameters.mes['links_names']:
            ret, lx = et.ENgetlinkindex(ln)
            link_index.append(lx)
        self.link_index = link_index

        pump_index = []
        for ln in self.sw_parameters.pumps['names']:
            ret, lx = et.ENgetlinkindex(ln)
            pump_index.append(lx)
        self.pumps_index = pump_index

    def get_node_index(self):
        mes_node_index = []
        print(sw_par.mes['nodes_names'])
        for nn in sw_par.mes['nodes_names']:
            ret, nx = et.ENgetnodeindex(nn)
            mes_node_index.append(nx)
            self.mes_node_index = mes_node_index
            print(nn, mes_node_index)
        aaaa

        tank_index = []
        for nn in self.sw_parameters.tanks['names']:
            ret, nx = et.ENgetnodeindex(nn)
            tank_index.append(nx)
            self.tank_index = tank_index

    def get_node_pattern(self):
        node_pattern_index= []
        for nn in sw_par.mes_pressure_number['inside']:
            print(self.mes_node_index[nn])
            ret, nx = et.ENgetnodevalue(self.mes_node_index[nn], et.EN_PATTERN)
            node_pattern_index.append(int(nx))
            print(node_pattern_index)
            self.node_pattern_index = node_pattern_index
        aaaaa
    def get_pumps_pattern_index(self):
        pattern_index = []
        for pn in self.sw_parameters.pumps['patterns_names']:
            ret, px = et.ENgetpatternindex(pn)
            pattern_index.append(px)
        self.pumps_pattern_index = pattern_index

    def get_node_pattern_names(self):
        node_pattern_names = []
        for pn in self.node_pattern_index:
            ret, px = et.ENgetpatternid(pn)
            node_pattern_names.append(px)
        self.node_pattern_names = node_pattern_names

    def get_pattern_values(self):
        pattern_value = []
        self.demand_patterns_values = []
        print(self.node_pattern_index)
        aaaaaa
        for pn in self.node_pattern_index:
            ret, T = et.ENgetpatternlen(pn)
            print(pn, T)
            for t in range(1,T+1):
                ret, px = et.ENgetpatternvalue(pn,t)
                pattern_value.append(px)
            self.demand_patterns_values.append(pattern_value)

    def set_time_duration(self):
        et.ENsettimeparam(0, self.sw_parameters.time['duration_s'])

    def set_tank_inital(self):
        for tank_in, initial_val in zip(self.tank_index, self.sw_parameters.tanks['initial_lev']):
            et.ENsetnodevalue(tank_in, 8, initial_val)

    def set_pattern_values(self,species, move_time):
        for i in range(self.sw_parameters.number['pumps']):
            for t in range(self.sw_parameters.number['hydraulic_steps']):
                 et.ENsetpatternvalue(self.pumps_pattern_index[i], t + 1, species[i*\
                                                                        self.sw_parameters.number['hydraulic_steps'] + t])

        for i in range(len(self.node_pattern_index)):
            for t in range(self.sw_parameters.number['hydraulic_steps']):
                et.ENsetpatternvalue(self.node_pattern_index[i], t + 1,
                                     self.demand_patterns_values[i][t + move_time])

    def save_temp_file(self):
        self.set_time_duration()

        et.ENsaveinpfile('temporary_' + sw_par.file_name[0])

    def get_set_parameters(self):
        self.open_epanet()
        self.get_link_index()
        self.get_node_index()
        self.get_pattern_index()
        self.set_time_duration()

        #self.save_temp_file()

    def prepare_empty_dict_to_comput(self):

        self.data = {'error_output': []}

        self.data['time_output'] = []

        for mes_node in self.sw_parameters.mes['nodes_names']:
            self.data['pressure_output_' + mes_node] = []

        for mes_link in self.sw_parameters.mes['links_names']:
            self.data['flow_output_' + mes_link] = []

        for mes_tank in self.sw_parameters.tanks['names']:
            self.data['tank_output_' + mes_tank] = []

        for pumps_name in self.sw_parameters.pumps['names']:
            self.data['energy_output_' + pumps_name] = []

        for tank_input in self.sw_parameters.tanks['names']:
            self.data['tank_input_' + tank_input] = []

        for pump_input in self.sw_parameters.pumps['names']:
            self.data['pump_input_' + pump_input] = []

        for demand_input in self.node_pattern_names:
            self.data['demand_input_' + demand_input] = []

    def get_hydraulic_values(self):

        et.ENopenH()
        et.ENinitH(0)
        time = []
        flow = [[] for i in range(self.sw_parameters.number['mes_links'])]
        pressure = [[] for i in range(len(self.node_index))]
        energy = [[] for i in range(self.sw_parameters.number['pumps'])]
        error = []
        while True:
            ret, t = et.ENrunH()

            if t % self.sw_parameters.time['hydraulic_step_s'] == 0:
                error.append(ret)
                time.append(t)

                if 'flow' in self.sw_parameters.hydraulic_monitor_values:
                    for i in range(0, len(self.link_index)):
                        ret, p = et.ENgetlinkvalue(self.link_index[i], et.EN_FLOW)
                        flow[i].append(p)

                if 'pressure' in self.sw_parameters.hydraulic_monitor_values:
                    for i in range(0, len(self.node_index)):
                        ret, p = et.ENgetnodevalue(self.node_index[i], et.EN_PRESSURE)
                        pressure[i].append(p)

                if 'energy' in self.sw_parameters.hydraulic_monitor_values:
                    for i in range(0, self.sw_parameters.number['pumps']):
                        ret, p = et.ENgetlinkvalue(self.pumps_index[i], et.EN_ENERGY)
                        energy[i].append(p)

            ret, tstep = et.ENnextH()
            if (tstep <= 0):
                break
        ret = et.ENcloseH()

        return flow, energy, pressure, error, time

    def insert_data(self, pressure, flow, energy, error, time):

        for ii in range(0, self.sw_parameters.number['mes_nodes']):
            self.data['pressure_output_' + self.sw_parameters.mes['nodes_names'][ii]].append(pressure[ii][:-1])

        for ii in range(self.sw_parameters.number['mes_nodes'] - self.sw_parameters.number['tanks'],
                        self.sw_parameters.number['mes_nodes']):
            self.data['tank_output_' + self.sw_parameters.tanks['names'][
                self.sw_parameters.number['mes_nodes'] - ii - 1]].append(pressure[ii][:-1])

        for ii in range(0, self.sw_parameters.number['tanks']):
            temp_tank_init = []
            for i in range(0, self.sw_parameters.time['duration_h']):
                temp_tank_init.append(self.sw_parameters.tanks['initial_lev'][ii])
            self.data['tank_input_' + self.sw_parameters.tanks['names'][ii]].append(temp_tank_init)

        for ii in range(self.sw_parameters.number['mes_links']):
            self.data['flow_output_' + self.sw_parameters.mes['links_names'][ii]].append(flow[ii][:-1])

        for ii in range(self.sw_parameters.number['pumps']):
            self.data['energy_output_' + self.sw_parameters.pumps['names'][ii]].append(energy[ii][:-1])

        self.data['time_output'].append(time[:-1])

        self.data['error_output'].append(error[:-1])

    def save_last_net(self,species):
        self.open_epanet()
        self.set_time_duration()
        self.set_tank_inital()
        self.set_pattern_values(species,0)
        self.save_temp_file()

    def get_data(self,species, move_time):

        self.prepare_empty_dict_to_comput()


        #self.open_epanet()
        self.set_time_duration()

        self.set_tank_inital()

        self.set_pattern_values(species, move_time)

        [flow, energy, pressure, error, time] = self.get_hydraulic_values()

        self.insert_data(pressure, flow, energy, error, time)

        self.close_epanet()

def main():
        print('')


if __name__ == '__main__':
    main()