"""

PCFLS_model.py

Run with python2.7+ (but not python3)

This file reads topology, traffic matrix, tunnel, and implements the FFC model in gurobi format.

"""

import sys
#sys.path.append('/opt/gurobi/new/lib/python2.7/')
from gurobipy import *
import numpy as np
from collections import defaultdict
import logging

def prepare_data(
    cap_file, tm_file, tm_index, unit_cap_file, tunnel_file, num_tunnel,
    is_symmetric, tunnel_node_disjoint, tunnel_edge_disjoint):

  # Index, Edge, Capacity definition
  nodes_set = set()
  arcs = []
  cap = {}
  with open(cap_file) as fin:
    for line in fin:
      if 'i' not in line:
        i, j, tcap = line.strip().split()
        tcap = float(tcap)
        nodes_set.add(i)
        nodes_set.add(j)
        if tcap > 0.0:
          arcs.append((i, j))
          cap[i,j] = tcap
  arcs_set = set(arcs)
  nodes = list(nodes_set)

  # Demand definition
  demand = {}
  sd_pairs = []
  with open(tm_file) as fin:
    for line in fin:
      if 's' not in line:
        s, t, h, tm = line.strip().split()
        tm = float(tm)
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
        i, j, u_cap = line.strip().split()
        u_cap = float(u_cap)
        if cap[i,j] % u_cap != 0:
          logging.warning(
              "cap[%s, %s] = %f is not integer multiple of ucap[%d, %d] = %f"
              % (i, j, cap[i,j], i, j, u_cap))
        n[i,j] = cap[i,j] // u_cap
        ucap[i,j] = u_cap

  # Tunnel definition
  tunnels = [str(t) for t in xrange(num_tunnel)]
  tunnels_set = set(tunnels)

  # Tunnel node and edge set definition
  tunnel_node_set, tunnel_edge_set = defaultdict(set), defaultdict(set)
  with open(tunnel_file) as fin:
    for line in fin:
      if 's' not in line:
        s, t, k, edge_list = line.strip().split()
        edge_list = edge_list.split(',')
        for e in edge_list:
          u, v = e.split('-')
          tunnel_edge_set[s,t,k].add((u, v))
          tunnel_node_set[s,t,k].add(u)
          tunnel_node_set[s,t,k].add(v)

  # Compute disjoint edge
  q = {}
  for s, t in sd_pairs:
    max_share = 0
    for i, j in arcs:
      g = 0
      for k in tunnels:
        if (i, j) in tunnel_edge_set[s,t,k] or (j, i) in tunnel_edge_set[s,t,k]:
          g = g + 1
      if g > max_share:
        max_share = g
    q[s, t] = max_share

  logging.info("q: %s" % q)
  logging.debug("tcap: %s" % cap)
  logging.debug("demand: %s" % demand)
  logging.debug("arcs: %s" % arcs)
  logging.debug("n: %s" % n)

  return {"nodes": nodes, "nodes_set": nodes_set, "arcs": arcs,
          "arcs_set": arcs_set, "tunnels": tunnels, "tunnels_set": tunnels_set,
          "capacity": cap, "unit_capacity": ucap, "demand": demand,
          "sd_pairs": sd_pairs, "sd_pairs_set": sd_pairs_set,
          "number_of_sublinks": n, "tunnel_node_set": tunnel_node_set,
          "tunnel_edge_set": tunnel_edge_set, "tunnel_node_disjoint":
          tunnel_node_disjoint, "tunnel_edge_disjoint": q}

def create_base_model(is_symmetric, data):
  # Parameters
  nodes, arcs, tunnels = data['nodes'], data['arcs'], data['tunnels']
  nodes_set, arcs_set, tunnels_set = \
      data['nodes_set'], data['arcs_set'], data['tunnels_set']
  sd_pairs = data['sd_pairs']
  cap, demand = data['capacity'], data['demand']
  tunnel_edge_set = data['tunnel_edge_set']

  # Gurobi Solver Model
  m = Model('tunnel_dual_method')

  # Variable definition
  mlu = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="mlu")
  a, b, pi, lambda_ = {}, {}, {}, {}
  for k in tunnels:
    for s, t in sd_pairs:
      a[k,s,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="a[%s,%s,%s]" % (k, s, t))
      pi[k,s,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="pi[%s,%s,%s]" % (k, s, t))

  for s, t in sd_pairs:
    b[s,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="b[%s,%s]" % (s, t))
    lambda_[s,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="lambda_[%s,%s]" % (s, t))

  # Objective definition
  m.setObjective(mlu, GRB.MINIMIZE)

  # Constraints definition
  for s, t in sd_pairs:
    m.addConstr(
        b[s,t] == demand[s,t], "c-2[%s,%s]" % (s, t)
    )

  for i, j in arcs:
    m.addConstr(
        quicksum(a[k,s,t] for k in tunnels for s, t in sd_pairs if (i, j) in
        tunnel_edge_set[s,t,k]) <= mlu * cap[i,j], "c-3[%s,%s]" % (i, j)
    )

  for k in tunnels:
    for s, t in sd_pairs:
      m.addConstr(
          pi[k,s,t] + lambda_[s,t] >= a[k,s,t], "dual-1[%s,%s,%s]" % (k, s, t)
      )

  m.update()
  return m

def compute_mlu(
    base_model, num_node_failure, num_link_failure, data, output_file):
  m = base_model.copy()

  # Data
  nodes, arcs, tunnels = data['nodes'], data['arcs'], data['tunnels']
  sd_pairs = data['sd_pairs']
  p, q = data['tunnel_node_disjoint'], data['tunnel_edge_disjoint']

  # Variable retrieval
  a, b, pi, lambda_ = {}, {}, {}, {}
  for k in tunnels:
    for s, t in sd_pairs:
      a[k,s,t] = m.getVarByName("a[%s,%s,%s]" % (k, s, t))
      pi[k,s,t] = m.getVarByName("pi[%s,%s,%s]" % (k, s, t))

  for s, t in sd_pairs:
    b[s,t] = m.getVarByName("b[%s,%s]" % (s, t))
    lambda_[s,t] = m.getVarByName("lambda_[%s,%s]" % (s, t))

  # Constraint definition
  for s, t in sd_pairs:
    m.addConstr(
        quicksum(a[k,s,t] for k in tunnels)
        - (num_node_failure * p + num_link_failure * q[s,t]) * lambda_[s,t]
        - quicksum(pi[k,s,t] for k in tunnels)
        >= b[s,t], "c-1[%s,%s]" % (s, t)
    )

  # Solve
  #m.Params.nodemethod = 2
  m.Params.method = 2
  #m.Params.Crossover = 0
  m.Params.Threads = 8
  m.optimize()

  logging.debug("Runtime: %f seconds" % m.Runtime)


  out_file = open(output_file, 'w')
 
  out_file.write('Physical tunnel reservation:\n')
  out_file.write('s t k r\n')
  for s, t in sd_pairs:
    for k in tunnels:
      out_file.write('%s %s %s %s\n' % (s, t, k, a[k,s,t].X))

  if m.Status == GRB.OPTIMAL or m.Status == GRB.SUBOPTIMAL:
    return m.ObjVal, m.Runtime
  return None, None
