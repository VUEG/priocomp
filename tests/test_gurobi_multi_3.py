#!/usr/bin/python

# Copyright 2016, Gurobi Optimization, Inc.

# Want to cover three different sets but subject to a common budget of
# elements allowed to be used. However, the sets have different priorities to
# be covered; and we tackle this by using multi-objective optimization.

from __future__ import print_function
from gurobipy import *


try:
    # Sample data
    Groundset = range(20)
    Fraction  = 0.2
    rep  = [3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 6, 0, 0, 0, 0, 0, 1, 0]
    cost = [2, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1]
    frac = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
    SetObjWeight = [1.0, -1.0]

    # Create initial model
    model = Model('ilp')

    # Initialize decision variables for ground set:
    vars = model.addVars(Groundset, vtype=GRB.BINARY, name='rEl')

    # Fraction constraint
    frac_constr = LinExpr()
    frac_constr.addTerms(frac, vars.values())
    # Match frac constraint EQUAL
    model.addConstr(lhs=frac_constr, sense=GRB.EQUAL, rhs=Fraction)

    # Set global sense for ALL objectives
    model.ModelSense = GRB.MAXIMIZE

    # Limit how many solutions to collect
    model.setParam(GRB.Param.PoolSolutions, 100)

    # Set number of objectives
    model.NumObj = 2

    # Objective 1: maximize representation
    rep_obj = LinExpr()
    rep_obj.addTerms(rep, vars.values())
    model.setParam(GRB.Param.ObjNumber, 0)
    model.ObjNWeight = 1.0
    model.ObjNName = 'Representation'
    model.ObjNRelTol = 0.01
    model.ObjNAbsTol = 1.0
    model.setAttr(GRB.Attr.ObjN, rep_obj)

    # Objective 2: minimize cost
    cost_obj = LinExpr()
    cost_obj.addTerms(cost, vars.values())
    model.setParam(GRB.Param.ObjNumber, 1)
    model.ObjNWeight = -1.0
    model.ObjNName = 'Cost'
    model.ObjNRelTol = 0.01
    model.ObjNAbsTol = 2.0
    model.setAttr(GRB.Attr.ObjN, cost_obj)

    model.update()
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
        print('Optimization was stopped with status ' + str(status))
        sys.exit(1)

    print("\nrep:    {}".format(rep))
    print("cost:   {}".format(cost))

     # Print best selected set
    print('Selected elements in best solution:')
    for e in Groundset:
        if vars[e].X > 0.9:
            print(' El%d' % e, end='')
    print('')

    # Print number of solutions stored
    nSolutions = model.SolCount
    print('Number of solutions found: ' + str(nSolutions))

    # Print objective values of solutions
    if nSolutions > 10:
        nSolutions = 10
    print('Objective values for first ' + str(nSolutions) + ' solutions:')
    for i in [0, 1]:
        model.setParam(GRB.Param.ObjNumber, i)
        print('\tSet%d' % i, end='')
        for e in range(nSolutions):
            model.setParam(GRB.Param.SolutionNumber, e)
            print(' %6g' % model.ObjNVal, end='')
        print('')

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))

except AttributeError as e:
    print('Encountered an attribute error: ' + str(e))
