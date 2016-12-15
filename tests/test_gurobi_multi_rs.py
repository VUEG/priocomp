#!/usr/bin/python
import numpy as np
import numpy.ma as ma
from gurobipy import *

# Helper functions ------------------------------------------------------------

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx]


# Set seed
np.random.seed(3)
# Use a common NoData value
nodata_value = -1

# Generate data ---------------------------------------------------------------
nrow = 100
ncol = 100
# Mask for NoData
mask = np.random.randint(0, 2, size=(nrow, ncol)).astype(np.bool)
# Representation matrix (to be maximized). Implemented as a masked array
rep = ma.masked_array(np.round(np.random.rand(nrow, ncol), 3), mask=mask)
# Cost matrix (to be minimized).
cost = ma.masked_array(np.round(np.random.rand(nrow, ncol) * 100, 1), mask=mask)
# Change the masked values to nodata_value
ma.set_fill_value(rep, nodata_value)
ma.set_fill_value(cost, nodata_value)

# Optimization ----------------------------------------------------------------

# Get all the non-masked data as a 1D arrays
rep_array = ma.compressed(rep)
cost_array = ma.compressed(cost)

# Fractions of elements out of the total number of elements
elem_fractions = np.full(rep_array.size, 1 / rep_array.size, dtype=np.float)
# Calculate cumsum to define what is the closest fraction multiple to
# desired fraction
cs_elem_fractions = np.cumsum(elem_fractions)
# Total fraction ([0, 1]) of the elements in 'values' that can be used.
# E.g. fraction of 0.2 corresponds to 4 elements when there are 20
# elements in total
fraction = 0.2
# Find the index of closest cumsum value to the requested fraction
fraction_approx = find_nearest(cs_elem_fractions, fraction)

try:

    # Create initial model
    model = Model('multiobj')

    # Initialize decision variables for ground set:
    for element in np.nditer(rep_array):
        model.addVar(vtype=GRB.BINARY)
    # Integrate new variables
    model.update()
    vars = model.getVars()

    # Constraint: limit total number of elements to be picked to be at most
    # fraction of total number of elements
    frac_constr = LinExpr()
    frac_constr.addTerms(elem_fractions.tolist(), vars)
    model.addConstr(lhs=frac_constr, sense=GRB.EQUAL, rhs=fraction_approx)

    # Set global sense for ALL objectives
    model.ModelSense = GRB.MAXIMIZE

    # Limit how many solutions to collect
    model.setParam(GRB.Param.PoolSolutions, 100)

    # Set number of objectives
    model.NumObj = 2

    # Objective 1: maximize values
    model.setParam(GRB.Param.ObjNumber, 0)
    model.ObjNWeight = 1
    model.ObjNName = "Values"
    #model.ObjNRelTol = 0.01
    #model.ObjNAbsTol = 1.0
    model.setAttr(GRB.Attr.ObjN, vars, rep_array.tolist())

    # Objective 2: minimize values
    model.setParam(GRB.Param.ObjNumber, 1)
    # Minimizing an objective function is equivalent to maximizing the
    # negation of that function -> use ObjNWeight = -1.0
    model.ObjNWeight = -1.0
    model.ObjNName = "Costs"
    #model.ObjNRelTol = 0.01
    #model.ObjNAbsTol = 2.0
    model.setAttr(GRB.Attr.ObjN, vars, cost_array.tolist())

    # Save problem
    # model.write('multiobj.lp')

    # Optimize
    model.optimize()

    model.setParam(GRB.Param.OutputFlag, 0)

    # Status checking
    status = model.Status
    if status == GRB.Status.INF_OR_UNBD or \
       status == GRB.Status.INFEASIBLE or \
       status == GRB.Status.UNBOUNDED:
        print('The model cannot be solved because it is infeasible or unbounded')
        print('Optimization was stopped with status {}'.format(status))
        sys.exit(1)

    if status != GRB.Status.OPTIMAL:
        print('Optimization was stopped with status {}'.format(status))
        sys.exit(1)

    # Print number of solutions stored
    nSolutions = model.SolCount
    print('\nNumber of solutions found: {}'.format(nSolutions))

    # Print objective values of solutions
    if nSolutions > 10:
        nSolutions = 10
    print('Objective values for first {} solutions:'.format(nSolutions))
    for i in [0, 1]:
        model.setParam(GRB.Param.ObjNumber, i)

        print('\t{0: >6}'.format(model.ObjNName), end='')
        for e in range(nSolutions):
            model.setParam(GRB.Param.SolutionNumber, e)
            print(' {0: >#014.2f}'.format(model.ObjNVal), end='')
        print('')

    res = np.asarray([var.x for var in model.getVars()], dtype=np.uint8)

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))

except AttributeError as e:
    print('Encountered an attribute error: ' + str(e))
