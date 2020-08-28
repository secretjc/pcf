#!/usr/bin/env python

"""

tunnel.py

Run with python2.7+ (but not python3)
"""

####
#### Imports
####
import sys
import logging

import module.pcf_models.tunnel_ffc_mlu as tunnel_ffc_mlu
import module.pcf_models.tunnel_dual_method_mlu as tunnel_dual_mlu
import module.pcf_models.tunnel_virtual_ab as tunnel_virtual_ab
import module.pcf_models.tunnel_virtual_ab_com as tunnel_virtual_ab_com

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

class TunnelFailureFFCMluSolver(object):
  def __init__(self, _sentinel=None, main_config=None, topo_config=None,
               solver_config=None):
    self.main_config = main_config
    self.topo_config = topo_config
    self.solver_config = solver_config
    self.logger = logging.getLogger("TunnelFailureFFCSolver")
    self.base_model = None
    self._parse_configs()

  def _parse_configs(self):
    assert self.topo_config is not None, "No topo config found!"
    self.cap_file = self.topo_config['data']['cap_file']
    self.unit_cap_file = self.topo_config['data']['cap_file']
    self.tm_file = self.topo_config['data']['tm_file']
    self.tunnel_file = self.topo_config['data']['tunnel_file']
    self.tunnel_node_disjoint = 0
    self.tunnel_edge_disjoint = 0
    self.num_tunnel = self.topo_config['attributes']['num_parallel_tunnels']
    self.B = 0
    self.tm_index = self.topo_config['traffic_matrix']['tm_index']
    self.output_path = self.main_config['output']
    #TODO
    self.is_symmetric = True

  def compute_max_throughput(self, _sentinel=None, num_link_failure=None):
    if self.base_model == None:
      self.logger.debug("No base model. Creating base model")
      self.data = tunnel_ffc_mlu.prepare_data(
          self.cap_file, self.tm_file, self.tm_index, self.unit_cap_file,
          self.tunnel_file, self.num_tunnel, self.is_symmetric,
          self.tunnel_node_disjoint, self.tunnel_edge_disjoint)
      self.base_model = tunnel_ffc_mlu.create_base_model(self.is_symmetric, self.data)

    return tunnel_ffc_mlu.compute_max_throughput(
        self.base_model, 0, num_link_failure, self.data, self.output_path)

class TunnelFailureMluSolver(object):
  def __init__(self, _sentinel=None, main_config=None, topo_config=None,
               solver_config=None):
    self.main_config = main_config
    self.topo_config = topo_config
    self.solver_config = solver_config
    self.logger = logging.getLogger("TunnelFailureSolver")
    self.base_model = None
    self._parse_configs()

  def _parse_configs(self):
    assert self.topo_config is not None, "No topo config found!"
    self.cap_file = self.topo_config['data']['cap_file']
    self.unit_cap_file = self.topo_config['data']['cap_file']
    self.tm_file = self.topo_config['data']['tm_file']
    self.tunnel_file = self.topo_config['data']['tunnel_file']
    self.num_tunnel = self.topo_config['attributes']['num_parallel_tunnels']
    self.B = 0
    self.tm_index = self.topo_config['traffic_matrix']['tm_index']
    self.output_path = self.main_config['output']
    #TODO
    self.is_symmetric = True

  def compute_max_throughput(self, _sentinel=None, num_link_failure=None):
    if self.base_model == None:
      self.logger.debug("No base model. Creating base model")
      self.data = tunnel_dual_mlu.prepare_data(
          self.cap_file, self.tm_file, self.tm_index, self.unit_cap_file,
          self.tunnel_file, self.num_tunnel, self.is_symmetric)
      self.base_model = tunnel_dual_mlu.create_base_model(self.is_symmetric, self.data)

    return tunnel_dual_mlu.compute_max_throughput(
        self.base_model, 0, num_link_failure, self.data, self.output_path)

class TunnelVirtualABSolver(object):
  def __init__(self, _sentinel=None, main_config=None, topo_config=None,
               solver_config=None):
    self.main_config = main_config
    self.topo_config = topo_config
    self.logger = logging.getLogger("TunnelVirtualNetworkSolver")
    self.base_model = None
    self._parse_configs()

  def _parse_configs(self):
    assert self.topo_config is not None, "No topo config found!"
    self.cap_file = self.topo_config['data']['cap_file']
    self.unit_cap_file = self.topo_config['data']['cap_file']
    self.tm_file = self.topo_config['data']['tm_file']
    self.tunnel_file = self.topo_config['data']['tunnel_file']
    self.B = 0
    self.tm_index = self.topo_config['traffic_matrix']['tm_index']
    self.output_path = self.main_config['output']
    self.is_edge_tunnel = False
    #TODO
    self.is_symmetric = True

  def compute_max_throughput(self, _sentinel=None, num_link_failure=None):
    if self.base_model == None:
      self.logger.debug("No base model. Creating base model")
      self.data = tunnel_virtual_ab.prepare_data(
          self.cap_file, self.tm_file, self.tm_index, self.unit_cap_file,
          self.tunnel_file, self.is_symmetric, self.is_edge_tunnel)
      self.base_model = tunnel_virtual_ab.create_base_model(self.is_symmetric, self.data)

    return tunnel_virtual_ab.compute_max_throughput(
        self.base_model, num_link_failure * 2, self.data, self.output_path)

class TunnelVirtualABComSolver(object):
  def __init__(self, _sentinel=None, main_config=None, topo_config=None,
               solver_config=None):
    self.main_config = main_config
    self.topo_config = topo_config
    self.logger = logging.getLogger("TunnelVirtualNetworkSolver")
    self.base_model = None
    self._parse_configs()

  def _parse_configs(self):
    assert self.topo_config is not None, "No topo config found!"
    self.cap_file = self.topo_config['data']['cap_file']
    self.unit_cap_file = self.topo_config['data']['cap_file']
    self.tm_file = self.topo_config['data']['tm_file']
    self.tunnel_file = self.topo_config['data']['tunnel_file']
    self.tm_index = self.topo_config['traffic_matrix']['tm_index']
    self.output_path = self.main_config['output']
    self.B = 0
    self.is_edge_tunnel = False
    #TODO
    self.is_symmetric = True

  def compute_max_throughput(self, _sentinel=None, num_link_failure=None):
    if self.base_model == None:
      self.logger.debug("No base model. Creating base model")
      self.data = tunnel_virtual_ab_com.prepare_data(
          self.cap_file, self.tm_file, self.tm_index, self.unit_cap_file,
          self.tunnel_file, self.is_symmetric, self.is_edge_tunnel)
      self.base_model = tunnel_virtual_ab_com.create_base_model(self.is_symmetric, self.data)

    return tunnel_virtual_ab_com.compute_max_throughput(
        self.base_model, num_link_failure * 2, self.data, self.output_path)
