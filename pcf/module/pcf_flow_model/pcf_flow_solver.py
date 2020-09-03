#!/usr/bin/env python

"""

pcf_flow_solver.py

Run with python2.7+ (but not python3)
"""

####
#### Imports
####
import sys
import logging

sys.path.append('../..')
import module.pcf_flow_model.restricted_PCF_flow_model as restricted_PCF_flow_model

####
#### Authorship information
####
__author__ = "Yiyang Chang"
__copyright__ = "Copyright 2016, Purdue ISL Robustness Framerwork Project"
__credits__ = ["Yiyang Chang"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Yiyang Chang"
__email__ = "chang256@purdue.edu"

class RestrictedPCFFlowSolver(object):
  def __init__(self, _sentinel=None, main_config=None, topo_config=None,
               solver_config=None):
    self.main_config = main_config
    self.topo_config = topo_config
    self.solver_config = solver_config
    self.logger = logging.getLogger("ModifiedR3Solver")
    self._parse_configs()
    self.base_model = None
    self.base_lb_model = None
    self.model_set = None
    self.u1 = None
    self.u2 = {}

  def _parse_configs(self):
    assert self.topo_config is not None, "No topo config found!"
    self.cap_file = self.topo_config['data']['cap_file']
    self.unit_cap_file = self.topo_config['data']['cap_file']
    self.tm_file = self.topo_config['data']['tm_file']
    self.B = 0
    self.tm_index = self.topo_config['traffic_matrix']['tm_index']
    self.f = int(self.topo_config['attributes']['num_link_failures'])

  def _fix_vars(self, prefix, is_symmetric, sorted_edges):
    if prefix is "":
      return []

    fix_edge_list = []

    for i, bit in enumerate(prefix.split('|')):
      if bit != '?':
        fix_edge_list.append((str(sorted_edges[i][0]), str(sorted_edges[i][1]), int(bit)))
        if is_symmetric:
          fix_edge_list.append((str(sorted_edges[i][1]), str(sorted_edges[i][0]), int(bit)))

    return fix_edge_list

  def compute_mlu(self, _sentinel=None, method=None,
                  is_symmetric=False, do_profile=False, write_lp=False,
                  search_stack=[], prefix='', use_epsilon=False,
                  verbose=False, sorted_edges=None, output_file="output.txt"):
    """
    """
    if self.base_model == None:
      self.logger.debug("No base model. Creating base model")
      self.base_model, self.model_set = \
          restricted_PCF_flow_model.create_base(self.cap_file, self.tm_file, self.tm_index,
                                  self.unit_cap_file, is_symmetric)
    self.logger.debug("Add additional constraints and solve")
    fix_edge_list = self._fix_vars(prefix, is_symmetric, sorted_edges)
    return restricted_PCF_flow_model.compute_mlu(self.base_model, self.f,
        fix_edge_list, self.model_set, is_symmetric, sorted_edges,
        set(search_stack), output_file)
