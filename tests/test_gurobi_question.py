#!/usr/bin/python

from gurobipy import *

# Values to be maximized
values = [3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 6, 0, 0, 0, 0, 0, 1, 0]
# Cost of each element in 'values'
cost =   [2, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1]
# Total budget for selecting items from 'values'
budget = 3
# Fractions of elements out of the total number of elements
elem_fractions = [1 / len(values)] * len(values)
# Total fraction ([0, 1]) of the elements in 'values' that can be used.
# E.g. fraction of 0.2 corresponds to 4 elements when there are 20
# elements in total
fraction = 0.2

# Create the model object
model = Model('ilp')

# Initialize decision variables:
vars = model.addVars(range(len(values)), vtype=GRB.BINARY, name='el')
model.update()

# Set objective: maximize values
obj = LinExpr()
obj.addTerms(values, vars.values())
model.setObjective(obj, sense=GRB.MAXIMIZE)

# Cost constraint
cost_constr = LinExpr()
cost_constr.addTerms(cost, vars.values())
model.addConstr(lhs=cost_constr, sense=GRB.LESS_EQUAL, rhs=budget)

# Fraction constraint
frac_constr = LinExpr()
frac_constr.addTerms(elem_fractions, vars.values())
model.addConstr(lhs=frac_constr, sense=GRB.EQUAL, rhs=fraction)

model.update()
# Optimize
model.optimize()

print("\nvalues: {}".format(values))
print("cost:   {}".format(cost))

# Print best selected set
print('\nSelected elements in best solution:')
obj_sum = 0
cost_sum = 0
for e in range(len(values)):
    if vars[e].X > 0.9:
        obj_sum += values[e]
        cost_sum += cost[e]
        print(' element {} (rep: {}, cost: {})'.format(e, values[e], cost[e]), end='\n')
print('')

print("Budget:    {}".format(budget))
print("Fraction:  {}".format(fraction))
print("Objective: {}".format(obj_sum))
print("Cost:      {}".format(cost_sum))