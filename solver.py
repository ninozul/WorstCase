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
import os, sys, shutil, csv

# GOComp modules - this should be visible on the GOComp evaluation system
import data

# modules for this code
#sys.path.append(os.path.normpath('.')) # better way to make this visible?
#import something

class Solver():

    def __init__(self):
        self.data = data.Data()

    def write_sol1_bus_section(self, w):

        w.writerow(['--bus section'])
        w.writerow(['i', 'v', 'theta', 'b'])
        for b in self.data.raw.buses.values():
            vhi = min(b.nvhi, b.evhi)
            vlo = max(b.nvlo, b.evlo)
            vm = 0.5 * (vhi + vlo)
            w.writerow([b.i, vm, 0.0, 0.0])

    def write_sol1_generator_section(self, w):

        w.writerow(['--generator section'])
        w.writerow(['i', 'uid', 'p', 'q'])
        for g in self.data.raw.generators.values():
            pg = 0.5 * (g.pt + g.pb) if g.stat == 1 else 0.0
            qg = 0.5 * (g.qt + g.qb) if g.stat == 1 else 0.0
            w.writerow([g.i, g.id, pg, qg])

    def write_sol1(self, sol_name):

        with open(sol_name, 'wb') as sol_file:
            w = csv.writer(sol_file, delimiter=",", quotechar="'", quoting=csv.QUOTE_MINIMAL)
            self.write_sol1_bus_section(w)
            self.write_sol1_generator_section(w)

    def write_sol2_bus_section(self, w, k):

        w.writerow(['--bus section'])
        w.writerow(['i', 'v', 'theta', 'b'])
        for b in self.data.raw.buses.values():
            vhi = min(b.nvhi, b.evhi)
            vlo = max(b.nvlo, b.evlo)
            vm = 0.5 * (vhi + vlo)
            w.writerow([b.i, vm, 0.0, 0.0])

    def write_sol2_generator_section(self, w, k):

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

    def write_sol2_delta_section(self, w, k):

        w.writerow(['--delta section'])
        w.writerow(['delta'])
        w.writerow([0.0])

    def write_sol2_ctg(self, w, k):

        w.writerow(['--contingency'])
        w.writerow(['label'])
        w.writerow([k.label])
        self.write_sol2_bus_section(w, k)
        self.write_sol2_generator_section(w, k)
        self.write_sol2_delta_section(w, k)

    def write_sol2(self, sol_name):

        with open(sol_name, 'wb') as sol_file:
            w = csv.writer(sol_file, delimiter=",", quotechar="'", quoting=csv.QUOTE_MINIMAL)
            for k in self.data.con.contingencies.values():
                self.write_sol2_ctg(w, k)

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



