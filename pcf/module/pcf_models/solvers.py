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

import module.pcf_models.FFC_model as FFC_model
import module.pcf_models.PCFTF_model as PCFTF_model
import module.pcf_models.PCFLS_model as PCFLS_model
import module.pcf_models.PCFLS_plus_model as PCFLS_plus_model

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

class FFC_Solver(object):
  def __init__(self, _sentinel=None, main_config=None, topo_config=None,
               solver_config=None):
    self.main_config = main_config
    self.topo_config = topo_config
    self.solver_config = solver_config
    self.logger = logging.getLogger("FFC_Solver")
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

  def compute_mlu(self, _sentinel=None, num_link_failure=None):
    if self.base_model == None:
      self.logger.debug("No base model. Creating base model")
      self.data = FFC_model.prepare_data(
          self.cap_file, self.tm_file, self.tm_index, self.unit_cap_file,
          self.tunnel_file, self.num_tunnel, self.is_symmetric,
          self.tunnel_node_disjoint, self.tunnel_edge_disjoint)
      self.base_model = FFC_model.create_base_model(self.is_symmetric, self.data)

    return FFC_model.compute_mlu(
        self.base_model, 0, num_link_failure, self.data, self.output_path)

class PCFTF_Solver(object):
  def __init__(self, _sentinel=None, main_config=None, topo_config=None,
               solver_config=None):
    self.main_config = main_config
    self.topo_config = topo_config
    self.solver_config = solver_config
    self.logger = logging.getLogger("PCFTF_Solver")
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

  def compute_mlu(self, _sentinel=None, num_link_failure=None):
    if self.base_model == None:
      self.logger.debug("No base model. Creating base model")
      self.data = PCFTF_model.prepare_data(
          self.cap_file, self.tm_file, self.tm_index, self.unit_cap_file,
          self.tunnel_file, self.num_tunnel, self.is_symmetric)
      self.base_model = PCFTF_model.create_base_model(self.is_symmetric, self.data)

    return PCFTF_model.compute_mlu(
        self.base_model, 0, num_link_failure, self.data, self.output_path)

class PCFLS_Solver(object):
  def __init__(self, _sentinel=None, main_config=None, topo_config=None,
               solver_config=None):
    self.main_config = main_config
    self.topo_config = topo_config
    self.logger = logging.getLogger("PCFLS_Solver")
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

  def compute_mlu(self, _sentinel=None, num_link_failure=None):
    if self.base_model == None:
      self.logger.debug("No base model. Creating base model")
      self.data = PCFLS_model.prepare_data(
          self.cap_file, self.tm_file, self.tm_index, self.unit_cap_file,
          self.tunnel_file, self.is_symmetric, self.is_edge_tunnel)
      self.base_model = PCFLS_model.create_base_model(self.is_symmetric, self.data)

    return PCFLS_model.compute_mlu(
        self.base_model, num_link_failure * 2, self.data, self.output_path)

class PCFLS_plus_Solver(object):
  def __init__(self, _sentinel=None, main_config=None, topo_config=None,
               solver_config=None):
    self.main_config = main_config
    self.topo_config = topo_config
    self.logger = logging.getLogger("PCFLS_plus_Solver")
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

  def compute_mlu(self, _sentinel=None, num_link_failure=None):
    if self.base_model == None:
      self.logger.debug("No base model. Creating base model")
      self.data = PCFLS_plus_model.prepare_data(
          self.cap_file, self.tm_file, self.tm_index, self.unit_cap_file,
          self.tunnel_file, self.is_symmetric, self.is_edge_tunnel)
      self.base_model = PCFLS_plus_model.create_base_model(self.is_symmetric, self.data)

    return PCFLS_plus_model.compute_mlu(
        self.base_model, num_link_failure * 2, self.data, self.output_path)
