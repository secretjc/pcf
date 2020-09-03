#!/usr/bin/env python

"""

generate_tunnels.py

Run with python2.7+ (but not python3)
"""
####
#### Imports
####
import sys
import time
import logging
import yaml
import argparse
from collections import defaultdict
import itertools
import random
import networkx as nx
from itertools import islice

#import local gurobi path if needed.
sys.path.append('/package/gurobi/8.0.1/lib/python2.7/')
sys.path.append('..')
import module.pcf_flow_model.pcf_flow_solver as pcf_flow_solver

####
#### Authorship information
####
__author__ = "Chuan Jiang"
__copyright__ = "Copyright 2020, Purdue ISL PCF Project"
__credits__ = ["Chuan Jiang"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Chuan Jiang"
__email__ = "jiang486@purdue.edu"
__status__ = "beta"

def _parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument("--topo_config", default="../_config/b4.yaml",
                      help="topology config file")
  parser.add_argument("--tunnel_type", default="physical")
  parser.add_argument("--output_path", default="tunnel.tab")
  return parser.parse_args()

def _parse_configs(args):

  with open(args.topo_config, 'r') as f:
    topo_config = yaml.load(f)

  return {'topo_config': topo_config}

def bellman(weight, s, t, n_node):
    d = {}
    prev = {}
    for i in range(n_node):
        d[i] = 1000
    d[s] = 0
    flag = True
    while flag:
        flag = False
        for (v1, v2) in weight:
            if d[v1] < 1000 and d[v1] + weight[(v1, v2)] < d[v2]:
                flag = True
                prev[v2] = v1
                d[v2] = d[v1] + weight[(v1, v2)]
    tunnel = []
    if d[t] < 1000:
        now = t
        while now != s:
            tunnel = [(prev[now],now)] + tunnel
            now = prev[now]
    return tunnel

def bhandari(edge, k, s, t, graph, n_node):
    weight = {}
    for u, v in edge:
        weight[(u, v)] = 1
    tunnels = []
    subgraph = set()
    n_disjoint_p = 0
    for i in range(k): #range(k)
        tunnel = bellman(weight, s, t, n_node)
        if tunnel != []:
            n_disjoint_p += 1
            tunnels.append(tunnel)
            for (u, v) in tunnel:
                if weight[(u, v)] == 1:
                    weight[(u, v)] = 2000
                    weight[(v, u)] = -1
                    subgraph.add((u,v))
                else:
                    weight[(u, v)] = 1
                    weight[(v, u)] = 1
                    subgraph.remove((v,u))  
    res = []
    new_g = nx.Graph()
    for i in range(n_node):
        new_g.add_node(i)
    for u, v in subgraph:
        new_g.add_edge(u, v)
    while n_disjoint_p > 0:
        path = nx.shortest_path(new_g, s, t)
        res.append(path)
        for i in range(len(path) - 1):
            new_g.remove_edge(path[i], path[i+1])
        n_disjoint_p -= 1
    #if len(res) == 1:
    #    print "ERROR", s, t
    #else:
    #    if not distinct(res[0], res[1]):
    #        print "ERRO", s, t
    if len(res) < k:
        k_paths = list(islice(nx.shortest_simple_paths(graph, s, t), k))
        for path in k_paths:
            if path not in res and len(res) < k:
                res.append(path)
    return res

def _to_path_string(path):
    path_str = str(path[0]) + "-"
    for j in range(len(path) - 2):
        path_str = path_str + str(path[j+1]) + "," + str(path[j+1]) + "-"
    path_str = path_str + str(path[len(path)-1])
    return path_str

def _to_path_string_skip(path):
    new_path = []
    for j in range(len(path)):
        if j % 2 == 0:
            new_path.append(path[j])
    path_str = str(new_path[0]) + "-"
    for j in range(len(new_path) - 2):
        path_str = path_str + str(new_path[j+1]) + "," + str(new_path[j+1]) + "-"
    path_str = path_str + str(new_path[len(new_path)-1])
    return path_str

def generate_tunnel(capa_file, output_file, flag):
    fin = open(capa_file)
    capa = {}
    edge = {}
    edge_ = {}
    share = {}
    count = 0
    n_node = 0
    for line in fin.readlines():
        if 'i' in line:
            line.strip()
        else:
            info = line.strip().split(' ')
            v1 = int(info[0])
            v2 = int(info[1])
            if v1 > n_node:
                n_node = v1
            if v2 > n_node:
                n_node = v2
            capa[(v1,v2)] = int(float(info[2])) * 1000
            edge[(v1,v2)] = count
            share[(v1,v2)] = 0
            edge_[count] = (v1,v2)
            count = count + 1
    n_node = n_node + 1

    demand = {}
    demand_ = {}
    count = 0
    for v1 in range(n_node):
        for v2 in range(n_node):
            if v1 != v2:
                demand[(v1,v2)] = 1
                demand_[count] = (v1,v2)
                count = count + 1
    graph = nx.Graph()
    for v in range(n_node):
        graph.add_node(v)
    for (u,v) in edge:
        graph.add_edge(u,v)
    k = 3
    a_path = []
    b_path = []
    paths = {}
    out_file = open(output_file, 'w')
    out_file.write("s t k edges\n")
    count = 0
    for s, t in demand:
            if s != t:
                count = count + 1
                pick = bhandari(edge, k, s, t, graph, n_node)
                #b_path.append((s,t,0,_to_path_string_skip(list(nx.shortest_path(graph, s, t)))))
                b_path.append((s,t,0,_to_path_string(list(nx.shortest_path(graph, s, t)))))
                for i in range(len(pick)):
                    path = pick[i]
                    path_str = _to_path_string(path)
    	            out_file.write("%s %s %s %s\n" % (s, t, i, path_str))
                    skip_path_str = _to_path_string_skip(path)
                    if i % 2 == 0:
                        a_path.append((s,t,i / 2,skip_path_str))
    if flag == 1:
        for path in b_path:
            s, t, i, path_str = path
            out_file.write("%s %s %s %s no\n" % (s, t, i, path_str))

def widest_path(e1, e2, weight, n_node):
    prev = {}
    cost = {}
    for i in range(n_node):
        prev[i] = -1
        cost[i] = -1
    cost[e1] = 10000000
    finished = set()
    while e2 not in finished:
        m = -1
        for i in range(n_node):
            if i not in finished and cost[i] > m:
                m = cost[i]
                now_node = i
        if m == -1:
            print 'ERROR'
        if m == -1 or now_node == e2:
            break
        for i in range(n_node):
            if (now_node, i) in weight and min(weight[(now_node, i)], cost[now_node]) > cost[i]:
                cost[i] = min(weight[(now_node, i)], cost[now_node])
                prev[i] = now_node
        finished.add(now_node)
    now_node = e2
    bottle_neck = cost[e2]
    path = []
    while now_node != -1:
        path = [now_node] + path
        if prev[now_node] != -1:
            weight[(prev[now_node], now_node)] -= bottle_neck
        now_node = prev[now_node]
    return _to_path_string(path)

def _generate_from_pcf_flow(capa_file, r3_file, output_file):
    n_node = 0
    fin = open(capa_file)
    edge = {}
    for line in fin.readlines():
        if 'i' in line:
            line.strip()
        else:
            info = line.strip().split(' ')
            v1 = int(info[0])
            v2 = int(info[1])
            if v1 + 1 > n_node:
                n_node = v1 + 1
            if v2 + 1 > n_node:
                n_node = v2 + 1
            edge[(v1,v2)] = {}
    fin.close()

    sd = {}
    for s in range(n_node):
        for t in range(n_node):
            if s != t:
                sd[(s,t)] = {}

    fin = open(r3_file)
    for line in fin.readlines():
        if 'p:' in line:
            info = line.strip().split(' ')
            v1 = int(info[1])
            v2 = int(info[2])
            e1 = int(info[3])
            e2 = int(info[4])
            w = float(info[5])
            edge[(e1,e2)][(v1,v2)] = w
        if 'r:' in line:
            info = line.strip().split(' ')
            v1 = int(info[1])
            v2 = int(info[2])
            s = int(info[3])
            t = int(info[4])
            w = float(info[5])
            sd[(s,t)][(v1,v2)] = w
    fin.close()
    p_tunnel = {}
    for e1, e2 in edge:
        ls1 = widest_path(e1, e2, edge[(e1,e2)], n_node)
        p_tunnel[(e1,e2)] = [(ls1, str(e1) + '-' + str(e2))]
    out_file = open(output_file, 'w')
    for s, t in sd:
        ls = widest_path(s, t, sd[(s,t)], n_node)
        out_file.write("%s %s %s %s no\n" % (s, t, 0, ls))
        count = 1
        if (s,t) in p_tunnel:
            for ls, hint in p_tunnel[(s,t)]:
                out_file.write("%s %s %s %s %s\n" % (s, t, count, ls, hint))
                count += 1
    out_file.close()
  
def _compute_pcf_flow(topo_config, pcf_flow_output):
  logger = logging.getLogger('compute_pcf_flow')
  restricted_pcf_flow_solver = pcf_flow_solver.RestrictedPCFFlowSolver(
      main_config=None,
      topo_config=topo_config,
      solver_config=None)

  solver = restricted_pcf_flow_solver
  throughput, index, solving_time = solver.compute_mlu(output_file=pcf_flow_output)
  logger.info("PCF Flow mlu: %s" % (throughput))

def _main(args, configs):
  topo_config = configs['topo_config']
  logging.basicConfig(level='INFO')
  logger = logging.getLogger('main')

  if args.tunnel_type == 'physical':
    generate_tunnel(topo_config['data']['cap_file'], args.output_path, 0)
  if args.tunnel_type == 'LS':
    generate_tunnel(topo_config['data']['cap_file'], args.output_path, 1)
  if args.tunnel_type == 'CLS':
    generate_tunnel(topo_config['data']['cap_file'], args.output_path, 0)
    _compute_pcf_flow(topo_config, "pcf_flow_output.txt")
    _generate_from_pcf_flow(topo_config['data']['cap_file'], "pcf_flow_output.txt", args.output_path)
  logger.info("Done.")

if __name__ == "__main__":
  args = _parse_args()
  configs = _parse_configs(args)
  _main(args, configs)
