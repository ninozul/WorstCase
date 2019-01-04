'''
solver.py
main code for the worst case feasible solver
'''

# sol1 format
'''
--bus section
i, v, theta, b
1,1.05999994278,0.0,-5.0
2,1.04499995708,-0.066464908421,0.0
--generator section
i, uid, p, q
1,'1 ',175.294113159,1.8048504591
'''

# sol2 format
'''
--contingency
label
G_1_1
--bus section
i, v, theta, b
1,1.05999994278,0.0,-5.0
2,1.04499995708,-0.0601769089699,0.0
--generator section
i, uid, p, q
1,1 ,0.0,0.0
--delta section
delta
0.0
--contingency
label
L_1_2_BL
--bus section
i, v, theta, b
1,1.05999994278,0.0,-5.0
2,1.04499995708,-0.0601769089699,0.0
--generator section
i, uid, p, q
1,1 ,176.474899292,9.35432052612
--delta section
delta
0.0
'''

# built in imports
import os, sys, shutil, csv, time, StringIO

# GOComp modules - this should be visible on the GOComp evaluation system
import data

# modules for this code
#sys.path.append(os.path.normpath('.')) # better way to make this visible?
#import something

def csv2string(data):
    si = StringIO.StringIO()
    cw = csv.writer(si)
    cw.writerows(data)
    return si.getvalue()

class Solver():

    def __init__(self):
        self.data = data.Data()

    def form_sol1_bus_section_rows(self):

        start_time = time.time()
        self.sol1_bus_section_rows = [
            [b.i, 0.5 * (min(b.nvhi, b.evhi) + max(b.nvlo, b.evlo)), 0.0, 0.0]
            for b in self.data.raw.buses.values()]
        end_time = time.time()
        print('form sol1 bus section rows time: %f' % (end_time - start_time))

    def form_sol1_generator_section_rows(self):

        start_time = time.time()
        self.sol1_generator_section_rows = [
            [g.i, g.id, (0.5 * (g.pt + g.pb) if g.stat == 1 else 0.0), (0.5 * (g.qt + g.qb) if g.stat == 1 else 0.0)]
            for g in self.data.raw.generators.values()]
        end_time = time.time()
        print('form sol1 generator section rows time: %f' % (end_time - start_time))

    def form_sol2_bus_section_rows(self):
        '''Note this is independent of contingency and therefore
        need not be done more than once, at the start of the solve.
        This will save considerable time'''

        start_time = time.time()
        self.sol2_bus_section_rows = [
            [b.i, 0.5 * (min(b.nvhi, b.evhi) + max(b.nvlo, b.evlo)), 0.0, 0.0]
            for b in self.data.raw.buses.values()]
        self.sol2_bus_section_rows_string = csv2string(self.sol2_bus_section_rows)
        end_time = time.time()
        print('form sol2 bus section rows time: %f' % (end_time - start_time))

    def form_sol2_generator_section_rows_base(self):

        start_time = time.time()
        g_list = list(self.data.raw.generators.values())
        g_num_map = {g_list[i]:i for i in range(len(g_list))}
        self.gi = [g.i for g in g_list]
        self.gid = [g.id for g in g_list]
        self.gp0 = [(0.5 * (g.pt + g.pb) if g.stat == 1 else 0.0) for g in g_list]
        self.gq0 = [(0.5 * (g.qt + g.qb) if g.stat == 1 else 0.0) for g in g_list]
        self.sol2_generator_section_rows_base = [
            [self.gi[i], self.gid[i], self.gp0[i], self.gq0[i]]
            for i in range(len(g_list))]
        self.k_g_num_out_map = {
            k:[g_num_map[self.data.raw.generators[(e.i, e.id)]]
                for e in k.generator_out_events]
            for k in self.data.con.contingencies.values()}
        end_time = time.time()
        print('form sol2 generator section rows base time: %f' % (end_time - start_time))

    def restore_sol2_generator_section_rows_base(self, k):
        '''values for generators that are out of serrvice have been zeroed out.
        This function restores these values.'''

        #start_time = time.time()
        for i in self.k_g_num_out_map[k]:
            self.sol2_generator_section_rows_base[i][2:4] = [self.gp0[i], self.gq0[i]]
        #end_time = time.time()
        #print('restore sol2 generator section rows base time: %f' % (end_time - start_time))

    def form_sol2_generator_section_rows_ctg(self, k):
        '''Note this depends on contingency and therefore needs to be done
        once at the start of each contingency.
        This is somewhat inefficient and could probably be improved.
        The inefficiency is acceptable as there are typically far fewer
        generators than buses, so the bus sections dominate the execution time.'''

        # method 1
        '''
        start_time = time.time()
        g_off_id = [(e.i, e.id) for e in k.generator_out_events] # cheap
        print(time.time() - start_time)
        g_off = [ # third most
            self.data.raw.generators[e]
            for e in g_off_id]
        print(time.time() - start_time)
        g_on = [ # second most
            g for g in self.data.raw.generators.values()
            if (g.i, g.id) not in g_off_id]
        print(time.time() - start_time)
        self.sol2_generator_section_rows = ( # most expensive
            [[g.i, g.id, (0.5 * (g.pt + g.pb) if g.stat == 1 else 0.0), (0.5 * (g.qt + g.qb) if g.stat == 1 else 0.0)] for g in g_on] +
            [[g.i, g.id, 0.0, 0.0] for g in g_off])
        end_time = time.time()
        print('form sol2 generator section rows ctg time (method 1): %f' % (end_time - start_time))
        '''

        # method 2
        #'''
        #start_time = time.time()
        #self.sol2_generator_section_rows = copy.deepcopy(self.sol2_generator_section_rows_base) - deep copy is expensive, just restore the values
        self.sol2_generator_section_rows = self.sol2_generator_section_rows_base
        for i in self.k_g_num_out_map[k]:
            self.sol2_generator_section_rows[i][2:4] = [0.0, 0.0]
        #end_time = time.time()
        #print('form sol2 generator section rows ctg time (method 2): %f' % (end_time - start_time))
        #'''

    def write_sol1_bus_section(self, w):

        '''
        start_time = time.time()
        w.writerow(['--bus section'])
        w.writerow(['i', 'v', 'theta', 'b'])
        for b in self.data.raw.buses.values():
            vhi = min(b.nvhi, b.evhi)
            vlo = max(b.nvlo, b.evlo)
            vm = 0.5 * (vhi + vlo)
            w.writerow([b.i, vm, 0.0, 0.0])
        end_time = time.time()
        print('write sol1 bus section time: %f' % (end_time - start_time))
        '''

        start_time = time.time()
        w.writerow(['--bus section'])
        w.writerow(['i', 'v', 'theta', 'b'])
        w.writerows(self.sol1_bus_section_rows)
        end_time = time.time()
        print('write sol1 bus section time: %f' % (end_time - start_time))

    def write_sol1_generator_section(self, w):

        '''
        start_time = time.time()
        w.writerow(['--generator section'])
        w.writerow(['i', 'uid', 'p', 'q'])
        for g in self.data.raw.generators.values():
            pg = 0.5 * (g.pt + g.pb) if g.stat == 1 else 0.0
            qg = 0.5 * (g.qt + g.qb) if g.stat == 1 else 0.0
            w.writerow([g.i, g.id, pg, qg])
        end_time = time.time()
        print('write sol1 generator section time: %f' % (end_time - start_time))
        '''

        start_time = time.time()
        w.writerow(['--generator section'])
        w.writerow(['i', 'uid', 'p', 'q'])
        w.writerows(self.sol1_generator_section_rows)
        end_time = time.time()
        print('write sol1 generator section time: %f' % (end_time - start_time))

    def write_sol1(self, sol_name):

        self.form_sol1_bus_section_rows()
        self.form_sol1_generator_section_rows()
        with open(sol_name, 'wb') as sol_file:
            w = csv.writer(sol_file, delimiter=",", quotechar="'", quoting=csv.QUOTE_MINIMAL)
            self.write_sol1_bus_section(w)
            self.write_sol1_generator_section(w)

    def write_sol2_bus_section(self, sol_file, w, k):

        # method 1
        '''
        start_time = time.time()
        w.writerow(['--bus section'])
        w.writerow(['i', 'v', 'theta', 'b'])
        for b in self.data.raw.buses.values():
            vhi = min(b.nvhi, b.evhi)
            vlo = max(b.nvlo, b.evlo)
            vm = 0.5 * (vhi + vlo)
            w.writerow([b.i, vm, 0.0, 0.0])
        end_time = time.time()
        print("write sol2 bus section time (method 1): %f" % (end_time - start_time))
        '''

        # method 2
        '''
        start_time = time.time()
        w.writerow(['--bus section'])
        w.writerow(['i', 'v', 'theta', 'b'])
        w.writerows(self.sol2_bus_section_rows)
        end_time = time.time()
        print("write sol2 bus section time (method 2): %f" % (end_time - start_time))
        '''

        # method 3
        #start_time = time.time()
        w.writerow(['--bus section'])
        w.writerow(['i', 'v', 'theta', 'b'])
        sol_file.write(self.sol2_bus_section_rows_string)
        #end_time = time.time()
        #print("write sol2 bus section time (method 3): %f" % (end_time - start_time))

    def write_sol2_generator_section(self, w, k):

        # method 1
        '''
        start_time = time.time()
        w.writerow(['--generator section'])
        w.writerow(['i', 'uid', 'p', 'q'])
        for g in self.data.raw.generators.values():
            if (g.i, g.id) in [(e.i, e.id) for e in k.generator_out_events]:
                pg = 0.0
                qg = 0.0
            else:
                pg = 0.5 * (g.pt + g.pb) if g.stat == 1 else 0.0
                qg = 0.5 * (g.qt + g.qb) if g.stat == 1 else 0.0
            w.writerow([g.i, g.id, pg, qg])
        end_time = time.time()
        print("write sol2 generator section time (method 1): %f" % (end_time - start_time))
        '''

        # method 2 - requires calling form rows method at start of ctg
        #'''
        #start_time = time.time()
        w.writerow(['--generator section'])
        w.writerow(['i', 'uid', 'p', 'q'])
        w.writerows(self.sol2_generator_section_rows)
        #end_time = time.time()
        #print("write sol2 generator section time (method 2): %f" % (end_time - start_time))
        #'''

    def write_sol2_delta_section(self, w, k):

        w.writerow(['--delta section'])
        w.writerow(['delta'])
        w.writerow([0.0])

    def write_sol2_ctg(self, sol_file, w, k):

        #start_time = time.time()
        w.writerow(['--contingency'])
        w.writerow(['label'])
        w.writerow([k.label])
        self.write_sol2_bus_section(sol_file, w, k)
        self.write_sol2_generator_section(w, k)
        self.write_sol2_delta_section(w, k)
        #end_time = time.time()
        #print("write sol2 ctg time: %f" % (end_time - start_time))

    def write_sol2(self, sol_name):

        start_time = time.time()
        self.form_sol2_bus_section_rows()
        self.form_sol2_generator_section_rows_base()
        num_ctg = len(self.data.con.contingencies.values())
        num_ctg_done = 0
        with open(sol_name, 'wb') as sol_file:
            w = csv.writer(sol_file, delimiter=",", quotechar="'", quoting=csv.QUOTE_MINIMAL)
            start_time = time.time()
            for k in self.data.con.contingencies.values():
                self.form_sol2_generator_section_rows_ctg(k) # todo not needed with sol2 gen method 1
                self.write_sol2_ctg(sol_file, w, k)
                self.restore_sol2_generator_section_rows_base(k) # todo not needed with sol2 gen method 1
                num_ctg_done += 1
                print("ctg... total: %u, done: %u, time elapsed: %u" % (num_ctg, num_ctg_done, time.time() - start_time))
        end_time = time.time()
        print("write sol2 time: %f" % (end_time - start_time))

    def read_data(self, raw_name, rop_name, inl_name, con_name):

        print 'reading data files'
        print 'reading raw file: %s' % raw_name
        self.data.raw.read(os.path.normpath(raw_name))
        print 'reading rop file: %s' % rop_name
        self.data.rop.read(os.path.normpath(rop_name))
        print 'reading inl file: %s' % inl_name
        self.data.inl.read(os.path.normpath(inl_name))
        print 'reading con file: %s' % con_name
        self.data.con.read(os.path.normpath(con_name))
        print "buses: %u" % len(self.data.raw.buses)
        print "loads: %u" % len(self.data.raw.loads)
        print "fixed_shunts: %u" % len(self.data.raw.fixed_shunts)
        print "generators: %u" % len(self.data.raw.generators)
        print "nontransformer_branches: %u" % len(self.data.raw.nontransformer_branches)
        print "transformers: %u" % len(self.data.raw.transformers)
        print "areas: %u" % len(self.data.raw.areas)
        print "switched_shunts: %u" % len(self.data.raw.switched_shunts)
        print "generator inl records: %u" % len(self.data.inl.generator_inl_records)
        print "generator dispatch records: %u" % len(self.data.rop.generator_dispatch_records)
        print "active power dispatch records: %u" % len(self.data.rop.active_power_dispatch_records)
        print "piecewise linear cost functions: %u" % len(self.data.rop.piecewise_linear_cost_functions)
        print 'contingencies: %u' % len(self.data.con.contingencies)

'''
def write_sol1_bus_section(d, w):

    w.writerow(['--bus section'])
    for b in d.raw.buses.values():
        w.writerow([b.i, 0.5 * (min(b.nvhi, b.evhi) + max(b.nvlo, b.evlo)), 0.0, 0.0])

def write_sol1_generator_section(d, w):

    w.writerow(['--generator section'])
    for g in d.raw.generators.values():
        w.writerow([g.i, g.id, 0.5 * (g.pt + g.pb), 0.5 * (g.qt + g.qb)])

def write_sol1(d, sol1_name):

    with open(sol1_name, 'wb') as sol_file:
        w = csv.writer(sol_file, delimiter=",", quotechar="'", quoting=csv.QUOTE_MINIMAL)
        write_sol1_bus_section(d, w)
        write_sol1_generator_section(d, w)

def read_data(raw_name, rop_name, inl_name, con_name):

    print 'reading data files'
    p = data.Data()
    print 'reading raw file: %s' % raw_name
    if raw_name is not None:
        p.raw.read(os.path.normpath(raw_name))
    print 'reading rop file: %s' % rop_name
    if rop_name is not None:
        p.rop.read(os.path.normpath(rop_name))
    print 'reading inl file: %s' % inl_name
    if inl_name is not None:
        p.inl.read(os.path.normpath(inl_name))
    print 'reading con file: %s' % con_name
    if con_name is not None:
        p.con.read(os.path.normpath(con_name))
    print "buses: %u" % len(p.raw.buses)
    print "loads: %u" % len(p.raw.loads)
    print "fixed_shunts: %u" % len(p.raw.fixed_shunts)
    print "generators: %u" % len(p.raw.generators)
    print "nontransformer_branches: %u" % len(p.raw.nontransformer_branches)
    print "transformers: %u" % len(p.raw.transformers)
    print "areas: %u" % len(p.raw.areas)
    print "switched_shunts: %u" % len(p.raw.switched_shunts)
    print "generator inl records: %u" % len(p.inl.generator_inl_records)
    print "generator dispatch records: %u" % len(p.rop.generator_dispatch_records)
    print "active power dispatch records: %u" % len(p.rop.active_power_dispatch_records)
    print "piecewise linear cost functions: %u" % len(p.rop.piecewise_linear_cost_functions)
    print 'contingencies: %u' % len(p.con.contingencies)
    return p
'''





'''
--contingency
label
G_1_1
--bus section
i, v, theta, b
1,1.05999994278,0.0,-5.0
2,1.04499995708,-0.0601769089699,0.0
--generator section
i, uid, p, q
1,1 ,0.0,0.0
--delta section
delta
0.0
--contingency
label
L_1_2_BL
--bus section
i, v, theta, b
1,1.05999994278,0.0,-5.0
2,1.04499995708,-0.0601769089699,0.0
--generator section
i, uid, p, q
1,1 ,176.474899292,9.35432052612
--delta section
delta
0.0
'''



