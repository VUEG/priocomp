#!/usr/bin/python

# Copyright 2016, Gurobi Optimization, Inc.

# Want to cover three different sets but subject to a common budget of
# elements allowed to be used. However, the sets have different priorities to
# be covered; and we tackle this by using multi-objective optimization.

from __future__ import print_function
from gurobipy import *
import tabulate


def normalize(x):
    return [float(i)/max(x) for i in x]


try:
    # Values to be maximized
    values = [1, 14, 1, 1, 1, 1, 1, 1, 1, 14, 0, 0, 0, 14, 0, 0, 4, 0, 0, 14]
    values = [i / 100 for i in values]
    values = normalize(values)
    # Cost of each element in 'values'
    cost =   [0, 2, 1, 0, 1, 1, 1, 1, 1, 2, 0, 1, 1, 2, 0, 1, 3, 1, 1, 2]
    cost = normalize(cost)
    # Fractions of elements out of the total number of elements
    elem_fractions = [1 / len(values)] * len(values)
    # Total fraction ([0, 1]) of the elements in 'values' that can be used.
    # E.g. fraction of 0.2 corresponds to 4 elements when there are 20
    # elements in total
    fraction = 0.2

    # Create initial model
    model = Model('multiobj')

    # Initialize decision variables for ground set:
    vars = model.addVars(range(len(values)), vtype=GRB.BINARY, name='El')

    # Constraint: limit total number of elements to be picked to be at most
    # fraction of total number of elements
    frac_constr = LinExpr()
    frac_constr.addTerms(elem_fractions, vars.values())
    model.addConstr(lhs=frac_constr, sense=GRB.EQUAL, rhs=fraction)

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
    model.ObjNRelTol = 0.01
    model.ObjNAbsTol = 1.0
    model.setAttr(GRB.Attr.ObjN, vars, values)

    # Objective 2: minimize values
    model.setParam(GRB.Param.ObjNumber, 1)
    # Minimizing an objective function is equivalent to maximizing the
    # negation of that function -> use ObjNWeight = -1.0
    model.ObjNWeight = -1.0
    model.ObjNName = "Costs"
    model.ObjNRelTol = 0.01
    model.ObjNAbsTol = 2.0
    model.setAttr(GRB.Attr.ObjN, vars, cost)

    # Save problem
    model.write('multiobj.lp')

    # Optimize
    model.optimize()

    model.setParam(GRB.Param.OutputFlag, 0)

    # Status checking
    status = model.Status
    if status == GRB.Status.INF_OR_UNBD or \
       status == GRB.Status.INFEASIBLE  or \
       status == GRB.Status.UNBOUNDED:
        print('The model cannot be solved because it is infeasible or unbounded')
        sys.exit(1)

    if status != GRB.Status.OPTIMAL:
        print('Optimization was stopped with status {}'.format(status))
        sys.exit(1)

    # Find out which elements are in the solution
    selected = []
    for e in range(len(values)):
        if vars[e].X > 0.9:
            selected.append("*")
        else:
            selected.append(" ")

    # Print best selected set
    data = []
    for i in range(len(values)):
        data.append([i, values[i], cost[i], selected[i]])
    header = ["index", "values", "cost", "selected"]
    print(tabulate.tabulate(data, headers=header, tablefmt="pipe",
                            floatfmt=".2f"))

    # Print number of solutions stored
    nSolutions = model.SolCount
    print('\nNumber of solutions found: {}'.format(nSolutions))

    # Print objective values of solutions
    if nSolutions > 10:
        nSolutions = 10
    print('Objective values for first {} solutions:'.format(nSolutions))
    for i in [0, 1]:
        model.setParam(GRB.Param.ObjNumber, i)

        for e in range(nSolutions):
            model.setParam(GRB.Param.SolutionNumber, e)
            print('\t{0: >6}'.format(model.ObjNName), end='')
            print(' {0: >6}'.format(model.ObjNVal), end='')
        print('')

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))

except AttributeError as e:
    print('Encountered an attribute error: ' + str(e))
