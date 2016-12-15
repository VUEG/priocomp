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
    Subsets   = range(2)
    Budget    = 3
    Set = [ [ 1, 4, 1, 1, 1, 1, 1, 1, 1, 4, 0, 0, 0, 4, 0, 0, 4, 0, 0, 4 ],
            [ 0, 2, 0, 0, 0, 1, 1, 1, 1, 2, 0, 0, 0, 2, 0, 1, 3, 1, 1, 2 ]]
    elem_fractions = [1 / len(Set[0])] * len(Set[0])
    fraction = 0.2
    SetObjWeight   = [1.0, -1.0]
    SetObjName = ["Value", "Cost"]

    # Create initial model
    model = Model('multiobj')

    # Initialize decision variables for ground set:
    # x[e] == 1 if element e is chosen for the covering.
    Elem = model.addVars(Groundset, vtype=GRB.BINARY, name='El')

    # Constraint: limit total number of elements to be picked to be at most
    # Budget
    # Fraction constraint
    frac_constr = LinExpr()
    frac_constr.addTerms(elem_fractions, Elem.values())
    # Match frac constraint EQUAL
    model.addConstr(lhs=frac_constr, sense=GRB.EQUAL, rhs=fraction)

    # Set global sense for ALL objectives
    model.ModelSense = GRB.MAXIMIZE

    # Limit how many solutions to collect
    model.setParam(GRB.Param.PoolSolutions, 100)

    # Set number of objectives
    model.NumObj = 4

    # Set and configure i-th objective
    for i in Subsets:
        model.setParam(GRB.Param.ObjNumber, i)
        model.ObjNWeight   = SetObjWeight[i]

        model.ObjNName = SetObjName[i]
        model.ObjNRelTol = 0.01
        model.ObjNAbsTol = 1.0 + i
        model.setAttr(GRB.Attr.ObjN, Elem, Set[i])

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
        print('Optimization was stopped with status ' + str(status))
        sys.exit(1)

    # Print best selected set
    print('Selected elements in best solution:')
    selected = []
    for e in Groundset:
        if Elem[e].X > 0.9:
            print(' El%d' % e, end='')
            selected.append("*")
        else:
            selected.append(" ")
    print('')

    i_str = " ".join(["{0: <2}".format(i) for i in range(20)])
    print("\nindex:    {}".format(i_str))
    s_str = " ".join(["{0: <2}".format(i) for i in selected])
    print("selected: {}".format(s_str))
    v_str = " ".join(["{0: <2}".format(i) for i in Set[0]])
    print("values:   {}".format(v_str))
    c_str = " ".join(["{0: <2}".format(i) for i in Set[1]])
    print("cost:     {}".format(c_str))

    # Print number of solutions stored
    nSolutions = model.SolCount
    print('\nNumber of solutions found: ' + str(nSolutions))

    # Print objective values of solutions
    if nSolutions > 10:
        nSolutions = 10
    print('Objective values for first ' + str(nSolutions) + ' solutions:')
    for i in Subsets:
        model.setParam(GRB.Param.ObjNumber, i)

        for e in range(nSolutions):
            model.setParam(GRB.Param.SolutionNumber, e)
            print('\t{}'.format(model.ObjNName), end='')
            print(' %6g' % model.ObjNVal, end='')
        print('')

except GurobiError as e:
    print('Error code ' + str(e.errno) + ": " + str(e))

except AttributeError as e:
    print('Encountered an attribute error: ' + str(e))
