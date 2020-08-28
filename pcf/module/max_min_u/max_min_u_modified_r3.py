import sys
sys.path.append('/opt/gurobi/new/lib/python2.7/')
from gurobipy import *
import numpy as np
from collections import defaultdict
import logging
import copy
import random
import heapq
import time
import json

def _get_index_random(m, sorted_edges, search_stack_set):
  s1 = set(range(len(sorted_edges)))
  #s2 = set(search_stack)
  random.seed(0)
  try:
    rnd_idx = random.choice(list(s1 - search_stack_set))
    logging.debug("Randomly picked an index: %s" % rnd_idx)
    return rnd_idx
  except IndexError:
    logging.warning("IndexError exception catched.")
    return -1

def _get_index_of_max_x(m, sorted_edges, search_stack_set, arcs, unfixed_arcs,
                        cap):
  metric = defaultdict(float)
  omlu = m.getVarByName("omlu").X
  gap = 0.1
  #gap = 1.0
  for i, j in arcs:
    mlu = m.getVarByName("mlu[%s,%s]" % (i, j)).X
    logging.debug("U[%s,%s] = %f" % (i, j, mlu))
    if mlu >= omlu * (1 - gap): # Filter high U_ij
      for i1, j1 in unfixed_arcs:
        measure = m.getVarByName("p[%s,%s,%s,%s]" % (i, j, i1, j1)).X * cap[i1,j1] / cap[i,j]
        logging.debug("p_l(e) * c_l / c_e [%s,%s,%s,%s]: %f" % (i, j, i1, j1, measure))
        metric[i1,j1] += measure
  logging.debug("metric: %s" % [(i1, j1, metric[i1, j1]) for i1, j1 in sorted(metric.keys())])

  max_ = -1
  index = -1
  i_ret, j_ret = -1, -1
  for i, edge in enumerate(sorted_edges):
    u, v = str(edge[0]), str(edge[1])
    indicator = metric[u, v]
    if max_ < indicator and i not in search_stack_set:
      max_ = indicator
      index = i
      i_ret, j_ret = u, v
  if max_ < 1e-6:
    logging.warning("Maximum metric is < 1e-6!")
  logging.debug("Picked index: %s, i.e., Link (%s, %s)" % (index, i_ret, j_ret))
  return index


def _get_index_of_max_x_old(m, sorted_edges, search_stack_set, arcs,
                            unfixed_arcs, _shadow):
  #dual_sol_d_e = defaultdict(float)
  dual_sol_d_l = defaultdict(float)

  #for i, j in unfixed_arcs:
  for i, j in arcs:
    sum_dual_sol = 0
    #denominator = m.getConstrByName("pi-1[%s,%s]" % (i, j)).Pi
    for i1, j1 in unfixed_arcs:
      dual_sol = m.getConstrByName("pi-2[%s,%s,%s,%s]" % (i, j, i1, j1)).Pi
      #logging.debug("dual value of [%s,%s,%s,%s]: %f" % (i, j, i1, j1, dual_sol))
      sum_dual_sol += dual_sol
      dual_sol_d_l[i1,j1] += dual_sol
    #dual_sol_d_e[i,j] += sum_dual_sol
  #logging.debug("sum_l x_{e,l} forall e: %s" % dual_sol_d_e)
  logging.debug("sum_e x_{e,l} forall l: %s" % dual_sol_d_l)

  MAX = -1
  index = -1
  for i, edge in enumerate(sorted_edges):
    u, v = str(edge[0]), str(edge[1])
    x = dual_sol_d_l[u, v]
    if x > MAX and i not in search_stack_set:
      MAX = x
      index = i
  logging.debug("Picked index: %s" % index)
  return index

def _get_index_of_max_x_old_2(m, sorted_edges, search_stack_set, arcs,
                              unfixed_arcs, _shadow):
  #dual_sol_d_e = defaultdict(float)
  dual_sol_d_l = defaultdict(float)

  #for i, j in unfixed_arcs:
  for i, j in arcs:
    sum_dual_sol = 0
    denominator = m.getConstrByName("pi-1[%s,%s]" % (i, j)).Pi
    logging.debug("denominator for (%s, %s): %s" % (i, j, denominator))
    for i1, j1 in unfixed_arcs:
      dual_sol = m.getConstrByName("pi-2[%s,%s,%s,%s]" % (i, j, i1, j1)).Pi
      logging.debug("dual value of [%s,%s,%s,%s]: %f" % (i, j, i1, j1, dual_sol))
      dual_sol /= denominator
      sum_dual_sol += dual_sol
      dual_sol_d_l[i1,j1] += dual_sol
    #dual_sol_d_e[i,j] += sum_dual_sol
  #logging.debug("sum_l x_{e,l} forall e: %s" % dual_sol_d_e)
  logging.debug("sum_e x_{e,l} forall l: %s" % dual_sol_d_l)

  MAX = -1
  index = -1
  for i, edge in enumerate(sorted_edges):
    u, v = str(edge[0]), str(edge[1])
    x = dual_sol_d_l[u, v]
    if x > MAX and i not in search_stack_set:
      MAX = x
      index = i
  logging.debug("Picked index: %s" % index)
  return index

def prepare_data(cap_file, tm_file, tm_index, unit_cap_file, is_symmetric):
  # Index, Edge, Capacity definition
  node_i = []
  arcs = []
  cap = {}
  with open(cap_file) as fin:
    for line in fin:
      if 'i' not in line:
        i = line.strip().split()[0]
        j = line.strip().split()[1]
        tcap = float(line.strip().split()[2])
        if i not in node_i:
          node_i.append(i)
        if j not in node_i:
          node_i.append(j)
        if tcap > 0.0:
          arcs.append((i,j))
          cap[i,j] = tcap
  arcs_set = set(arcs)

  # Demand definition
  demand = defaultdict(float)
  sd_pairs = []
  with open(tm_file) as fin:
    for line in fin:
      if 's' not in line:
        s = line.strip().split()[0]
        t = line.strip().split()[1]
        h = line.strip().split()[2]
        tm = float(line.strip().split()[3])
        if int(h) == tm_index:
          demand[s,t] = tm
          if tm > 0.0:
            sd_pairs.append((s,t))
  sd_pairs_set = set(sd_pairs)

  # Sublink definition
  n = {}
  ucap = {}
  with open(unit_cap_file) as fin:
    for line in fin:
      if 'i' not in line:
        i = line.strip().split()[0]
        j = line.strip().split()[1]
        u_cap = float(line.strip().split()[2])
        if cap[i,j] % u_cap != 0:
          logging.warning(
              "cap[%s, %s] = %f is not integer multiple of ucap[%d, %d] = %f"
              % (i, j, cap[i,j], i, j, u_cap))
        n[i,j] = cap[i,j] // u_cap
        ucap[i,j] = u_cap

  logging.debug("node_i: %s" % node_i)
  logging.debug("tcap: %s" % cap)
  logging.debug("demand: %s" % demand)
  logging.debug("arcs: %s" % arcs)
  logging.debug("n: %s" % n)

  return [node_i, arcs, arcs_set, cap, ucap, demand, n, sd_pairs, sd_pairs_set]

def compute_mlu(base_model, f, fix_edge_list, model_set, is_symmetric,
                sorted_edges, search_stack_set, output_file):
  m = base_model.copy()

  if is_symmetric:
    f *= 2
  logging.debug("Number of failures : %d" % f)
  logging.debug("fix_edge_list: %s" % fix_edge_list)

  node_s, node_t, arcs, arcs_set, cap, ucap, demand, n, sd_pairs = model_set

  # Variable retrieval
  r = {}
  pi = {}
  ld = {}
  p = {}
  a = {}
  mlu = {}
  for i, j in arcs:
    for s, t in sd_pairs:
      r[i,j,s,t] = m.getVarByName("r[%s,%s,%s,%s]" % (i, j, s, t))

    for i1, j1 in arcs:
      pi[i,j,i1,j1] = m.getVarByName("pi[%s,%s,%s,%s]" % (i, j, i1, j1))
      p[i,j,i1,j1] = m.getVarByName("p[%s,%s,%s,%s]" % (i, j, i1, j1))

    a[i,j] = m.getVarByName("a[%s,%s]" % (i, j))
    ld[i,j] = m.getVarByName("ld[%s,%s]" % (i, j))
    mlu[i,j] = m.getVarByName("mlu[%s,%s]" % (i, j))

  beta = {}
  fixed_arcs_set = set([(i, j) for i,j,k in fix_edge_list])
  unfixed_arcs_set = arcs_set - fixed_arcs_set
  unfixed_arcs = list(unfixed_arcs_set)
  total_fixed = 0
  for i1, j1, k in fix_edge_list:
    total_fixed += k
    beta[i1,j1] = k
    for i, j in arcs:
      m.remove(m.getConstrByName("pi-2[%s,%s,%s,%s]" % (i, j, i1, j1)))

  #TODO: the p variable under this needs to be taken care of
  # constraint pi-1 sum((s,t), d[s,t] * r[i,j,s,t]) +
  #                 sum((i1,j1), pi[i,j,i1,j1]) + ld[i,j] * f <= mlu[i,j] *
  #                 c[i,j]
  # for i, j in arcs:
  #   m.addConstr(quicksum(demand[s,t] * r[i,j,s,t] for s, t in sd_pairs) +
  #       quicksum(n[i1,j1] * pi[i,j,i1,j1] for (i1, j1) in arcs_set) + ld[i,j] *
  #       f <= mlu[i,j] * cap[i,j], "pi-1[%s,%s]" % (i, j))
  for i, j in arcs:
    if (i, j) in unfixed_arcs_set:
      # m.addConstr(quicksum(demand[s,t] * r[i,j,s,t] for s, t in sd_pairs) +
      #     quicksum(n[i1,j1] * pi[i,j,i1,j1] for (i1, j1) in unfixed_arcs)
      #     + ld[i,j] * (f - total_fixed) <= mlu[i,j] * cap[i,j] - quicksum(k *
      #     p[i,j,i1,j1] * cap[i1,j1] for (i1, j1, k) in fix_edge_list),
      #     "pi-1[%s,%s]" % (i, j))
      m.addConstr(
          -quicksum(demand[s,t] * r[i,j,s,t] for s, t in sd_pairs)
              - quicksum(n[i1,j1] * pi[i,j,i1,j1] for (i1, j1) in unfixed_arcs)
              - ld[i,j] * (f - total_fixed)
              >=
              -mlu[i,j] * cap[i,j]
              + quicksum(k * p[i,j,i1,j1] for (i1, j1, k) in fix_edge_list),
          "pi-1[%s,%s]" % (i, j))
    else:
      # m.addConstr(quicksum(demand[s,t] * r[i,j,s,t] for s, t in sd_pairs) +
      #     quicksum(n[i1,j1] * pi[i,j,i1,j1] for (i1, j1) in unfixed_arcs)
      #     + ld[i,j] * (f - total_fixed) <= mlu[i,j] * ucap[i,j] * (n[i,j] - beta[i,j]) +
      #     a[i,j] * beta[i,j] - quicksum(k * p[i,j,i1,j1] * cap[i1,j1] for (i1,
      #     j1, k) in fix_edge_list), "pi-1[%s,%s]" % (i, j))
      m.addConstr(
          -quicksum(demand[s,t] * r[i,j,s,t] for s, t in sd_pairs)
              - quicksum(n[i1,j1] * pi[i,j,i1,j1] for (i1, j1) in unfixed_arcs)
              - ld[i,j] * (f - total_fixed)
              >=
              -mlu[i,j] * ucap[i,j] * (n[i,j] - beta[i,j]) - a[i,j] * beta[i,j]
              + quicksum(k * p[i,j,i1,j1] for (i1, j1, k) in fix_edge_list),
          "pi-1[%s,%s]" % (i, j))
    #else:
    #  m.addConstr(quicksum(demand[s,t] * r[i,j,s,t] for s, t in sd_pairs) <=
    #      mlu[i,j] * cap[i,j], "pi-1[%s,%s]" % (i, j))

  # Solve
  m.Params.method = 2
  m.Params.Crossover = 0
  m.Params.Threads = 6
  #m.Params.Presolve = 0
  #m.update()
  m.optimize()

  logging.info("Runtime: %f seconds" % m.Runtime)


  out_file = open(output_file, 'w')
  for i,j,i1,j1 in p:
    out_file.write("p: %s %s %s %s %s\n" % (i, j, i1, j1, p[i,j,i1,j1].X))

  for i,j,s,t in r:
    if r[i,j,s,t].X > 0: #demand[s,t] / 100:
      out_file.write("r: %s %s %s %s %s\n" % (i, j, s, t, r[i,j,s,t].X))

  if m.Status == GRB.OPTIMAL:
    if sorted_edges == None:
      return m.ObjVal, -1, m.Runtime
    else:
      return m.ObjVal, _get_index_of_max_x_old(m, sorted_edges,
        search_stack_set, arcs, unfixed_arcs, cap), m.Runtime
  else:
    logging.warning("Result is not optimal! Status = %s" % m.Status)
    if sorted_edges == None:
      return float('inf'), -1, m.Runtime
    else:
      return float('inf'), _get_index_random(m, sorted_edges, search_stack_set), m.Runtime

  return None

def create_base(cap_file, tm_file, tm_index, unit_cap_file, is_symmetric):
  # Index, Edge, Capacity definition
  node_i = []
  arcs = []
  cap = {}
  with open(cap_file) as fin:
    for line in fin:
      if 'i' not in line:
        i = line.strip().split()[0]
        j = line.strip().split()[1]
        tcap = float(line.strip().split()[2])
        if i not in node_i:
          node_i.append(i)
        if j not in node_i:
          node_i.append(j)
        if tcap > 0.0:
          arcs.append((i,j))
          cap[i,j] = tcap
  node_j = node_i
  node_s = node_i
  node_t = node_i
  arcs_set = set(arcs)

  demand = defaultdict(float)
  sd_pairs = []
  with open(tm_file) as fin:
    for line in fin:
      if 's' not in line:
        s = line.strip().split()[0]
        t = line.strip().split()[1]
        h = line.strip().split()[2]
        tm = float(line.strip().split()[3])
        if int(h) == tm_index:
          demand[s,t] = tm
          if tm > 0.0:
            sd_pairs.append((s,t))
  sd_pairs_set = set(sd_pairs)

  # Sublink definition
  n = {}
  ucap = {}
  with open(unit_cap_file) as fin:
    for line in fin:
      if 'i' not in line:
        i = line.strip().split()[0]
        j = line.strip().split()[1]
        u_cap = float(line.strip().split()[2])
        if cap[i,j] % u_cap != 0:
          logging.warning(
              "cap[%s, %s] = %f is not integer multiple of ucap[%d, %d] = %f"
              % (i, j, cap[i,j], i, j, u_cap))
        n[i,j] = cap[i,j] // u_cap
        ucap[i,j] = u_cap

  logging.debug("node_i: %s" % node_i)
  logging.debug("tcap: %s" % cap)
  logging.debug("demand: %s" % demand)
  logging.debug("arcs: %s" % arcs)
  logging.debug("n: %s" % n)

  m = Model('max_min_u_modified_r3')

  # Variable definition
  r = {}
  p = {}
  pi = {}
  ld = {}
  a = {}
  mlu = {}
  theta = {}
  for i, j in arcs:
    for s, t in sd_pairs:
      r[i,j,s,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS,
                            name="r[%s,%s,%s,%s]" % (i, j, s, t))

    for i1, j1 in arcs:
      p[i,j,i1,j1] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS,
                              name="p[%s,%s,%s,%s]" % (i, j, i1, j1))
      pi[i,j,i1,j1] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS,
                               name="pi[%s,%s,%s,%s]" % (i, j, i1, j1))
      if is_symmetric:
        theta[i,j,i1,j1] = m.addVar(vtype=GRB.CONTINUOUS,
                                    name="theta[%s,%s,%s,%s]" % (i, j, i1, j1))

    ld[i,j] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="ld[%s,%s]" % (i, j))
    a[i,j] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="a[%s,%s]" % (i, j))
    mlu[i,j] = m.addVar(vtype=GRB.CONTINUOUS, name="mlu[%s,%s]" % (i, j))

  omlu = m.addVar(vtype=GRB.CONTINUOUS, name="omlu")
  #sigma = m.addVar(lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="sigma")

  # Objective definition
  m.setObjective(omlu, GRB.MINIMIZE)
  #m.setObjective(omlu + 1 - sigma, GRB.MINIMIZE)

  # Constraints definition
  # constraint r-1 sum(j, r[s,j,s,t]) = 1
  for s, t in sd_pairs:
    m.addConstr(quicksum(r[s,j,s,t] for j in node_j if (s, j) in arcs_set)
        == 1, "r-1[%s,%s]" % (s, t))
    # m.addConstr(quicksum(r[s,j,s,t] for j in node_j if (s, j) in arcs_set)
    #     == sigma, "r-1[%s,%s]" % (s, t))

  # constraint r-2 r[j,s,s,t] = 1
  for j, s in arcs:
    for t in node_t:
      if (s, t) in sd_pairs_set:
        m.addConstr(r[j,s,s,t] == 0, "r-2[%s,%s,%s]" % (j, s, t))

  # constraint r-3 sum(j, r[i,j,s,t]) - sum(j, r[j,i,s,t]) = 0
  for i in node_i:
    for s, t in sd_pairs:
      if i != s and i != t:
        m.addConstr(quicksum(r[i,j,s,t] for j in node_j if (i, j) in
            arcs_set) - quicksum(r[j,i,s,t] for j in node_j if (j, i) in
            arcs_set) == 0, "r-3[%s,%s,%s]" % (i, s, t))

  # constraint p-1 sum(j, p[i1,j,i1,j1]) = a[i1,j]
  for i1, j1 in arcs:
    #TODO: This may have bug
    # m.addConstr(quicksum(p[i1,j,i1,j1] * cap[i1,j1] for j in node_j if (i1, j)
    #     in arcs_set) == a[i1,j1], "p-1[%s,%s]" % (i1, j1))
    # m.addConstr(quicksum(p[i1,j,i1,j1] * ucap[i1,j1] for j in node_j if (i1, j)
    #     in arcs_set) == a[i1,j1], "p-1[%s,%s]" % (i1, j1))
    m.addConstr(
        quicksum(p[i1,j,i1,j1] for j in node_j if (i1, j) in arcs_set)
            == a[i1,j1],
        "p-1[%s,%s]" % (i1, j1))

  # constraint p-2 p[j,i1,i1,j1] = 0
  for j, i1 in arcs:
    for j1 in node_j:
      if (i1, j1) in arcs:
        m.addConstr(p[j,i1,i1,j1] == 0, "p-2[%s,%s,%s]" % (j, i1, j1))

  # constraint p-3 sum(j, p[i,j,i1,j1]) - sum(j, p[j,i,i1,j1]) = 0
  for i1, j1 in arcs:
    for i in node_i:
      if i != i1 and i != j1:
        m.addConstr(quicksum(p[i,j,i1,j1] for j in node_j if (i, j) in
            arcs_set) - quicksum(p[j,i,i1,j1] for j in node_j if (j, i) in
            arcs_set) == 0, "p-3[%s,%s,%s]" % (i, i1, j1))

  # constraint pi-2 pi[i,j,i1,j1] + ld[i,j] >= p[i,j,i1,j1] +
  #                                            mlu[i,j] * (c[i,j] - a[i,j])
  for i, j in arcs:
    for i1, j1 in arcs:
      if is_symmetric:
        if i == i1 and j == j1:
          # The following works
          # m.addConstr(pi[i,j,i1,j1] + pi[i,j,j1,i1] + 2 * ld[i,j] >=
          #     p[i,j,i1,j1] * cap[i1,j1] + p[i,j,j1,i1] * cap[j1,i1] +
          #     mlu[i,j] * ucap[i,j] - a[i,j], "pi-2[%s,%s,%s,%s]" % (i, j, i1,
          #     j1))

          # The following also works
          #TODO: This may have bug
          # m.addConstr(pi[i,j,i1,j1] + ld[i,j] + theta[i,j,i1,j1] -
          #             theta[i,j,j1,i1] >= p[i,j,i1,j1] * cap[i1,j1] +
          #             mlu[i,j] * ucap[i,j] - a[i,j],
          #             "pi-2[%s,%s,%s,%s]" % (i, j, i1, j1))
          # m.addConstr(pi[i,j,i1,j1] + ld[i,j] + theta[i,j,i1,j1] -
          #             theta[i,j,j1,i1] >= p[i,j,i1,j1] * ucap[i1,j1] +
          #             mlu[i,j] * ucap[i,j] - a[i,j],
          #             "pi-2[%s,%s,%s,%s]" % (i, j, i1, j1))
          m.addConstr(
              pi[i,j,i1,j1] + ld[i,j] + theta[i,j,i1,j1] - theta[i,j,j1,i1]
                  >= p[i,j,i1,j1] + mlu[i,j] * ucap[i,j] - a[i,j],
              "pi-2[%s,%s,%s,%s]" % (i, j, i1, j1))

        else:
          # The following works
          # m.addConstr(pi[i,j,i1,j1] + pi[i,j,j1,i1] + 2 * ld[i,j] >=
          #     p[i,j,i1,j1] * cap[i1,j1] + p[i,j,j1,i1] * cap[j1,i1],
          #     "pi-2[%s,%s,%s,%s]" % (i, j, i1, j1))

          # The following also works
          #TODO: This may have bug
          # m.addConstr(pi[i,j,i1,j1] + ld[i,j] + theta[i,j,i1,j1] -
          #             theta[i,j,j1,i1] >= p[i,j,i1,j1] * cap[i1,j1],
          #             "pi-2[%s,%s,%s,%s]" % (i, j, i1, j1))
          # m.addConstr(pi[i,j,i1,j1] + ld[i,j] + theta[i,j,i1,j1] -
          #             theta[i,j,j1,i1] >= p[i,j,i1,j1] * ucap[i1,j1],
          #             "pi-2[%s,%s,%s,%s]" % (i, j, i1, j1))
          m.addConstr(
              pi[i,j,i1,j1] + ld[i,j] + theta[i,j,i1,j1] - theta[i,j,j1,i1]
                  >= p[i,j,i1,j1], "pi-2[%s,%s,%s,%s]" % (i, j, i1, j1))

      else:
        if i == i1 and j == j1:
          #TODO: This may have bug
          # m.addConstr(pi[i,j,i1,j1] + ld[i,j] >= p[i,j,i1,j1] * cap[i1,j1] +
          #     mlu[i,j] * ucap[i,j] - a[i,j], "pi-2[%s,%s,%s,%s]" % (i, j, i1,
          #     j1))
          m.addConstr(
              pi[i,j,i1,j1] + ld[i,j]
                  >= p[i,j,i1,j1] + mlu[i,j] * ucap[i,j] - a[i,j],
              "pi-2[%s,%s,%s,%s]" % (i, j, i1, j1))

        else:
          #TODO: This may have bug
          # m.addConstr(pi[i,j,i1,j1] + ld[i,j] >= p[i,j,i1,j1] * cap[i1,j1],
          #     "pi-2[%s,%s,%s,%s]" % (i, j, i1, j1))
          m.addConstr(
              pi[i,j,i1,j1] + ld[i,j] >= p[i,j,i1,j1],
              "pi-2[%s,%s,%s,%s]" % (i, j, i1, j1))

  # constraint mlu-1 omlu >= mlu[i,j]
  for i, j in arcs:
    m.addConstr(omlu >= mlu[i,j], "mlu-1[%s,%s]" % (i, j))

  # # constraint mlu-2 mlu[i,j] = mlu[j,i]
  # if is_symmetric:
  #   for i, j in arcs:
  #     m.addConstr(mlu[i,j] == mlu[j,i], "mlu-2[%s,%s]" % (i, j))

  m.update()
  return m, [node_s, node_t, arcs, arcs_set, cap, ucap, demand, n, sd_pairs]

