'''
MyPython1.py

Entry to Python worst case feasible solver, timed portion.
Syntax:

$ python MyPython1.py CON INL RAW ROP TIMELIMIT SCOREMETHOD NETWORKMODEL
'''

# built in imports
import sys, os

# modules for this code
sys.path.append(os.path.normpath('.')) # better way to make this visible?
import solver

args = sys.argv

con_name = args[1]
inl_name = args[2]
raw_name = args[3]
rop_name = args[4]
time_limit = args[5]
score_method = args[6]
network_model = args[7]
sol1_name = 'solution1.txt'
sol2_name = 'solution2.txt'

print('\nPython worst case feasible solver')
print('syntax:')
print('$ python MyPython1.py CON INL RAW ROP TIMELIMIT SCOREMETHOD NETWORKMODEL')
print('args:')
print(args)

# set up solver
s = solver.Solver()

# read the data files
s.read_data(raw_name, rop_name, inl_name, con_name)

# write the solution
s.write_sol1(sol1_name)
