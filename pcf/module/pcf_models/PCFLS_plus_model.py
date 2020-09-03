import sys
from gurobipy import *
import numpy as np
from collections import defaultdict
import logging

def prepare_data(
    cap_file, tm_file, tm_index, unit_cap_file, tunnel_file, is_symmetric, is_edge_tunnel):

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

  node_pairs = []
  atunnels = {}
  atraverse = {}
  btunnels = {}
  btraverse = {}
  impact = {}
  demand = {}
  for i in nodes_set:
    for j in nodes_set:
      if i != j:
        node_pairs.append((i, j))
      atunnels[i,j] = set()
      atraverse[i,j] = set()
      btunnels[i,j] = set()
      btraverse[i,j] = set()
      impact[i,j] = set()
      demand[i,j] = 0

  # Demand definition
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

  # Tunnel node and edge set definition
  atunnel_edge_set = {}
  btunnel_edge_set = {}
  btunnel_reason_set = {}
  logging.debug("tunnel_file: %s" % tunnel_file)
  anum_tunnel = 0
  bnum_tunnel = 0
  with open(tunnel_file) as fin:
    for line in fin:
      if 's' not in line:
        unpack = line.strip().split()
        if len(unpack) == 4:
	  #handle a tunnel
          s, t, k, edge_list = unpack
          edge_list = edge_list.split(',')
          atunnel_edge_set[anum_tunnel] = set()
          atunnels[s,t].add(anum_tunnel)
          for e in edge_list:
            u, v = e.split('-')
            atunnel_edge_set[anum_tunnel].add((u, v))
            atraverse[u,v].add(anum_tunnel)
          anum_tunnel = anum_tunnel + 1
        else:
          #handle b tunnel
          s, t, k, edge_list, reason_list = line.strip().split()
          edge_list = edge_list.split(',')
          btunnel_edge_set[bnum_tunnel] = set()
          btunnels[s,t].add(bnum_tunnel)
          for e in edge_list:
            u, v = e.split('-')
            btunnel_edge_set[bnum_tunnel].add((u, v))
            btraverse[u,v].add(bnum_tunnel)
          btunnel_reason_set[bnum_tunnel] = set()
          if "no" not in reason_list:
            reason_list = reason_list.split(',')
            for e in reason_list:
              u, v = e.split('-')
              btunnel_reason_set[bnum_tunnel].add((u, v))
              impact[u,v].add(bnum_tunnel)
          bnum_tunnel = bnum_tunnel + 1
  
  if is_edge_tunnel:
    for i,j in arcs:
      atunnels[i,j].add(anum_tunnel)
      atraverse[i,j].add(anum_tunnel)
      atunnel_edge_set[anum_tunnel] = set()
      atunnel_edge_set[anum_tunnel].add((i,j))
      anum_tunnel = anum_tunnel + 1
      btunnels[i,j].add(bnum_tunnel)
      btraverse[i,j].add(bnum_tunnel)
      impact[i,j].add(bnum_tunnel)
      btunnel_edge_set[bnum_tunnel] = set()
      btunnel_edge_set[bnum_tunnel].add((i,j))
      btunnel_reason_set[bnum_tunnel] = set()
      btunnel_reason_set[bnum_tunnel].add((i,j))
      bnum_tunnel = bnum_tunnel + 1

  logging.debug("tcap: %s" % cap)
  logging.debug("demand: %s" % demand)
  logging.debug("arcs: %s" % arcs)
  logging.debug("node_pairs: %s" % node_pairs)
  logging.debug("atunnel_edge_set: %s" % atunnel_edge_set)
  logging.debug("atunnels: %s" % atunnels)
  logging.debug("atraverse: %s" % atraverse)
  logging.debug("anum_tunnel: %s" % anum_tunnel)
  logging.debug("btunnel_edge_set: %s" % btunnel_edge_set)
  logging.debug("btunnel_reason_set: %s" % btunnel_reason_set)
  logging.debug("btunnels: %s" % btunnels)
  logging.debug("btraverse: %s" % btraverse)
  logging.debug("impact: %s" % impact)
  logging.debug("bnum_tunnel: %s" % bnum_tunnel)
  logging.debug("n: %s" % n)

  return {"nodes": nodes, "nodes_set": nodes_set, "arcs": arcs,
          "arcs_set": arcs_set, "capacity": cap, "unit_capacity": ucap, "demand": demand,
          "sd_pairs": sd_pairs, "node_pairs": node_pairs,
          "number_of_sublinks": n, 
          "atunnels": atunnels, "atraverse": atraverse, 
          "atunnel_edge_set": atunnel_edge_set, "anum_tunnel": anum_tunnel,
          "btunnels": btunnels, "btraverse": btraverse, "impact": impact, 
          "btunnel_edge_set": btunnel_edge_set, "btunnel_reason_set": btunnel_reason_set, 
          "bnum_tunnel": bnum_tunnel}

def create_base_model(is_symmetric, data):
  # Parameters
  nodes, arcs = data['nodes'], data['arcs']
  nodes_set, arcs_set = \
      data['nodes_set'], data['arcs_set']
  sd_pairs, node_pairs = data['sd_pairs'], data['node_pairs']
  cap, demand = data['capacity'], data['demand']
  atunnels = data['atunnels']
  atraverse = data['atraverse']
  atunnel_edge_set = data['atunnel_edge_set']
  anum_tunnel = data['anum_tunnel']
  btunnels = data['btunnels']
  btraverse = data['btraverse']
  impact = data['impact']
  btunnel_edge_set = data['btunnel_edge_set']
  btunnel_reason_set = data['btunnel_reason_set']
  bnum_tunnel = data['bnum_tunnel']

  # Gurobi Solver Model
  m = Model('tunnel_virtual_network')

  # Variable definition
  r, a, a0, b0, c0, p, pi, lambda_, amu, bmu, theta, asigma, bsigma, aphi, bphi = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
  for l in range(anum_tunnel):
    for t in nodes:
      r[l,t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="r[%s,%s]" % (l,t))
    a0[l] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="a0[%s]" % (l))
    for i, j in node_pairs:
      if l in atunnels[i,j]:
        amu[l,i,j] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="amu[%s,%s,%s]" % (l, i, j))
        aphi[l,i,j] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="aphi[%s,%s,%s]" % (l, i, j))
  for l in range(bnum_tunnel):
    b0[l] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="b0[%s]" % (l))
    c0[l] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="c0[%s]" % (l))
    for i, j in node_pairs:
      if l in btunnels[i,j] or l in btraverse[i,j]:
        bmu[l,i,j] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="bmu[%s,%s,%s]" % (l, i, j))
        bphi[l,i,j] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="bphi[%s,%s,%s]" % (l, i, j))
        for e1, e2 in btunnel_reason_set[l]:
          bsigma[l,i,j,e1,e2] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="bsigma[%s,%s,%s,%s,%s]" % (l, i, j, e1, e2))

  for i, j in node_pairs:
    lambda_[i,j] = m.addVar(vtype=GRB.CONTINUOUS, name="lambda_[%s,%s]" % (i, j))
    for e1, e2 in arcs:
      pi[i,j,e1,e2] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="pi[%s,%s,%s,%s]" % (i, j, e1, e2))
      theta[i,j,e1,e2] = m.addVar(vtype=GRB.CONTINUOUS, name="theta[%s,%s,%s,%s]" % (i, j, e1, e2))

  mlu = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name="mlu")

  # Objective definition
  m.setObjective(mlu, GRB.MINIMIZE)

  # Constraints definition
  for i, j in node_pairs:
    for e1, e2 in arcs:
      uc = 0
      m.addConstr(
        pi[i,j,e1,e2] + theta[i,j,e1,e2] - theta[i,j,e2,e1] + lambda_[i,j] + 
        quicksum((bsigma[l,i,j,e1,e2] - bphi[l,i,j]) for l in impact[e1,e2] if l in btunnels[i, j] or l in btraverse[i, j])+ 
        quicksum(- aphi[l,i,j] for l in atraverse[e1,e2] if l in atunnels[i, j]) >=
        uc,
        "dual_x[%s,%s,%s,%s]" % (i, j, e1, e2)
      )

  for i, j in node_pairs:
    uc = 0
    if (i,j) in arcs:
      uc = cap[i,j] * mlu
    m.addConstr(
      quicksum(a0[l] for l in atraverse[i,j]) <= uc,
      "tunnel_reserve[%s,%s]" % (i,j) 
    )

  for i, j in node_pairs:
    for l in atunnels[i,j]:
      m.addConstr(
        amu[l,i,j] + aphi[l,i,j] >= a0[l], "dual_ay[%s,%s,%s]" % (l,i,j)
      )

  for l in range(bnum_tunnel):
    for i, j in node_pairs:
      if l in btunnels[i,j] or l in btraverse[i,j]:
        rhs = 0
        if l in btraverse[i,j]:
          rhs = rhs - b0[l] + c0[l]
        if l in btunnels[i,j]:
          rhs = rhs + b0[l] - c0[l]
        m.addConstr(
          bmu[l,i,j] - quicksum(bsigma[l,i,j,e1,e2] for e1,e2 in btunnel_reason_set[l]) + bphi[l,i,j] >=
          rhs, "dual_by[%s,%s,%s]" % (l,i,j)
        )

  for l in range(anum_tunnel):
    m.addConstr(
      a0[l] >= 0, "cstr_a0[%s]" % (l)
    )

  for l in range(bnum_tunnel):
    m.addConstr(
      b0[l] >= 0, "cstr_b0[%s]" % (l)
    )

  for l in range(bnum_tunnel):
    m.addConstr(
      c0[l] >= 0, "cstr_c0[%s]" % (l)
    )

  m.update()
  return m

def compute_mlu(
    base_model, num_link_failure, data, output_file):
  m = base_model.copy()

  # Data
  nodes, arcs = data['nodes'], data['arcs']
  atunnels = data['atunnels']
  atraverse = data['atraverse']
  atunnel_edge_set = data['atunnel_edge_set']
  anum_tunnel = data['anum_tunnel']
  btunnels = data['btunnels']
  btraverse = data['btraverse']
  impact = data['impact']
  btunnel_edge_set = data['btunnel_edge_set']
  btunnel_reason_set = data['btunnel_reason_set']
  bnum_tunnel = data['bnum_tunnel']
  node_pairs = data['node_pairs']
  cap, demand = data['capacity'], data['demand']

  # Variable retrieval
  mlu = m.getVarByName("mlu")
  r, pi, lambda_, amu, bmu, a0, b0, c0, aphi, bphi, a, p, asigma, bsigma, theta = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}

  for l in range(anum_tunnel):
    a0[l] = m.getVarByName("a0[%s]" % (l))
    for i, j in node_pairs:
      amu[l,i,j] = m.getVarByName("amu[%s,%s,%s]" % (l, i, j))
      aphi[l,i,j] = m.getVarByName("aphi[%s,%s,%s]" % (l, i, j))
  for l in range(bnum_tunnel):
    b0[l] = m.getVarByName("b0[%s]" % (l))
    c0[l] = m.getVarByName("c0[%s]" % (l))
    for i, j in node_pairs:
      bmu[l,i,j] = m.getVarByName("bmu[%s,%s,%s]" % (l, i, j))
      bphi[l,i,j] = m.getVarByName("bphi[%s,%s,%s]" % (l, i, j))

  for i, j in node_pairs:
    lambda_[i,j] = m.getVarByName("lambda_[%s,%s]" % (i, j))
    for e1, e2 in arcs:
      pi[i,j,e1,e2] = m.getVarByName("pi[%s,%s,%s,%s]" % (i, j, e1, e2))
      theta[i,j,e1,e2] = m.getVarByName("theta[%s,%s,%s,%s]" % (i, j, e1, e2))

  # Constraint definition
  for i, j in node_pairs:
    m.addConstr(
      demand[i,j] + quicksum(pi[i,j,e1,e2] for e1,e2 in arcs) +
      lambda_[i,j] * num_link_failure + 
      quicksum(amu[l,i,j] for l in atunnels[i,j]) +
      quicksum(bmu[l,i,j] for l in range(bnum_tunnel) if l in btunnels[i,j] or l in btraverse[i,j]) -
      quicksum(a0[l] for l in atunnels[i,j]) -
      quicksum(b0[l] for l in btunnels[i,j]) + 
      quicksum(b0[l] for l in btraverse[i,j]) <= 0,
      "dual_obj[%s,%s]" % (i,j) 
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
  out_file.write('s t k r\n')
  for i, j in node_pairs:
    tunnel_list = list(atunnels[i,j])
    tunnel_list.sort()
    for k in range(len(tunnel_list)):
      out_file.write('%s %s %s %s\n' % (i, j, k, a0[tunnel_list[k]].X))

  out_file.write('Logical sequence reservation:\n')
  out_file.write('s t k r\n')
  for i, j in node_pairs:
    tunnel_list = list(btunnels[i,j])
    tunnel_list.sort()
    for k in range(len(tunnel_list)):
      out_file.write('%s %s %s %s\n' % (i, j, k, b0[tunnel_list[k]].X))

  out_file.write('Logical sequence(complement hint) reservation:\n')
  out_file.write('s t k r\n')
  for i, j in node_pairs:
    tunnel_list = list(btunnels[i,j])
    tunnel_list.sort()
    for k in range(len(tunnel_list)):
      out_file.write('%s %s %s %s\n' % (i, j, k, c0[tunnel_list[k]].X))

  if m.Status == GRB.OPTIMAL or m.Status == GRB.SUBOPTIMAL:
    return m.ObjVal, m.Runtime
  return None, None
