import sys
#sys.path.append('/opt/gurobi/new/lib/python2.7/')
from gurobipy import *
#import numpy as np
from collections import defaultdict
import logging

def prepare_data(
    cap_file, tm_file, tm_index, unit_cap_file, tunnel_file, num_tunnel,
    is_symmetric):

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

  print num_tunnel

  # Tunnel node and edge set definition
  tunnel_node_set, tunnel_edge_set = defaultdict(set), defaultdict(set)
  logging.debug("tunnel_file: %s" % tunnel_file)
  with open(tunnel_file) as fin:
    for line in fin:
      if 's' not in line:
        s, t, k, edge_list = line.strip().split()
        edge_list = edge_list.split(',')
        for e in edge_list:
          u, v = e.split('-')
          tunnel_edge_set[s,t,k].add((u, v))
          #tunnel_edge_set[s,t,k].add((v, u))
          tunnel_node_set[s,t,k].add(u)
          tunnel_node_set[s,t,k].add(v)


  logging.debug("tcap: %s" % cap)
  logging.debug("demand: %s" % demand)
  logging.debug("arcs: %s" % arcs)
  logging.debug("tunnel_node_set: %s" % tunnel_node_set)
  logging.debug("tunnel_edge_set: %s" % tunnel_edge_set)
  logging.debug("n: %s" % n)

  return {"nodes": nodes, "nodes_set": nodes_set, "arcs": arcs,
          "arcs_set": arcs_set, "tunnels": tunnels, "tunnels_set": tunnels_set,
          "capacity": cap, "unit_capacity": ucap, "demand": demand,
          "sd_pairs": sd_pairs, "sd_pairs_set": sd_pairs_set,
          "number_of_sublinks": n, "tunnel_node_set": tunnel_node_set,
          "tunnel_edge_set": tunnel_edge_set}

def create_base_model(is_symmetric, data):
  # Parameters
  nodes, arcs, tunnels = data['nodes'], data['arcs'], data['tunnels']
  nodes_set, arcs_set, tunnels_set = \
      data['nodes_set'], data['arcs_set'], data['tunnels_set']
  sd_pairs = data['sd_pairs']
  cap, demand = data['capacity'], data['demand']
  tunnel_node_set = data['tunnel_node_set']
  tunnel_edge_set = data['tunnel_edge_set']

  # Gurobi Solver Model
  m = Model('tunnel_dual_method')

  # Variable definition
  mlu = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="mlu")
  a, b, pi, lambda_, mu, theta, sigma, phi, beta = {}, {}, {}, {}, {}, {}, {}, {}, {}
  for k in tunnels:
    for s, t in sd_pairs:
      a[k,s,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="a[%s,%s,%s]" % (k, s, t))
      pi[k,s,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="pi[%s,%s,%s]" % (k, s, t))
      phi[k,s,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="phi[%s,%s,%s]" % (k, s, t))

  for s, t in sd_pairs:
    b[s,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="b[%s,%s]" % (s, t))
    mu[s,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="mu[%s,%s]" % (s, t))
    for i, j in arcs:
      sigma[i,j,s,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="sigma[%s,%s,%s,%s]" % (i, j, s, t))
      theta[i,j,s,t] = m.addVar(vtype=GRB.CONTINUOUS, name="theta[%s,%s,%s,%s]" % (i, j, s, t))


  # Objective definition
  #m.setObjective(quicksum(b[s,t] for s, t in sd_pairs), GRB.MAXIMIZE)
  m.setObjective(mlu, GRB.MINIMIZE)

  # Constraints definition
  #for s, t in sd_pairs:
  #  m.addConstr(
  #      b[s,t] <= demand[s,t], "c-2[%s,%s]" % (s, t)
  #  )

  for k in tunnels:
    for s, t in sd_pairs:
      if ((s,t,k) not in tunnel_edge_set):
        m.addConstr(a[k,s,t] == 0, "aempty[%s,%s,%s]" % (k, s, t))
        m.addConstr(pi[k,s,t] == 0, "piempty[%s,%s,%s]" % (k, s, t))
        m.addConstr(phi[k,s,t] == 0, "phiempty[%s,%s,%s]" % (k, s, t))

  for i, j in arcs:
    m.addConstr(
        quicksum(a[k,s,t] for k in tunnels for s, t in sd_pairs if (i, j) in
        tunnel_edge_set[s,t,k]) <= mlu * cap[i,j],
        "c-3[%s,%s]" % (i, j)
    )

  for k in tunnels:
    for s, t in sd_pairs:
      m.addConstr(
          pi[k,s,t] + phi[k,s,t] >= a[k,s,t], "dual-1[%s,%s,%s]" % (k, s, t)
      )

  for s, t in sd_pairs:
    for i, j in arcs:
      m.addConstr(
          -quicksum(pi[k,s,t] for k in tunnels if (i, j) in
          tunnel_edge_set[s,t,k]) + mu[s,t] + sigma[i,j,s,t] + theta[i,j,s,t] - theta[j,i,s,t] >= 0,
          "dual-3[%s,%s,%s,%s]" % (i, j, s, t)
      )

  m.update()
  return m

def compute_max_throughput(
    base_model, num_node_failure, num_link_failure, data, output_file):
  m = base_model.copy()

  # Data
  nodes, arcs, tunnels = data['nodes'], data['arcs'], data['tunnels']
  sd_pairs = data['sd_pairs']
  demand = data['demand']

  # Variable retrieval
  Z = m.getVarByName("mlu") 
  a, b, lambda_, mu, theta, sigma, phi = {}, {}, {}, {}, {}, {}, {}
  for k in tunnels:
    for s, t in sd_pairs:
      a[k,s,t] = m.getVarByName("a[%s,%s,%s]" % (k, s, t))
      phi[k,s,t] = m.getVarByName("phi[%s,%s,%s]" % (k, s, t))

  for s, t in sd_pairs:
    mu[s,t] = m.getVarByName("mu[%s,%s]" % (s, t))

  for i, j in arcs:
    for s, t in sd_pairs:
      sigma[i,j,s,t] = m.getVarByName("sigma[%s,%s,%s,%s]" % (i, j, s, t))
      theta[i,j,s,t] = m.getVarByName("theta[%s,%s,%s,%s]" % (i, j, s, t))

  # Constraint definition
  for s, t in sd_pairs:
    m.addConstr(
        quicksum(a[k,s,t] for k in tunnels)
        - (num_link_failure * 2 * mu[s,t]
           + quicksum(sigma[i,j,s,t] for i, j in arcs)
           + quicksum(phi[k,s,t] for k in tunnels))
        >= demand[s,t], "c-1[%s,%s]" % (s, t)
    )

  # Solve
  #m.Params.nodemethod = 2
  m.Params.method = 2
  #m.Params.Crossover = 0
  m.Params.Threads = 4
  m.optimize()

  logging.debug("Runtime: %f seconds" % m.Runtime)

  out_file = open(output_file, 'w') 
  out_file.write('Physical tunnel reservation:\n')
  for s, t in sd_pairs:
    for k in tunnels:
      out_file.write('%s %s %s %s\n' % (s, t, k, a[k,s,t].X))

  if m.Status == GRB.OPTIMAL or m.Status == GRB.SUBOPTIMAL:
    return m.ObjVal, m.Runtime
  return None, None
