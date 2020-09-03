#!/usr/bin/env python

"""

run.py

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
#import local gurobi path if needed.
sys.path.append('/package/gurobi/8.0.1/lib/python2.7/')
sys.path.append('..')
import module.pcf_models.solvers as solvers

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
  parser.add_argument("--main_config", default="../_config/main.yaml",
                      help="main config file")
  parser.add_argument("--topo_config", default="../_config/b4_config.yaml",
                      help="topology config file")
  return parser.parse_args()

def _parse_configs(args):
  with open(args.main_config, 'r') as f:
    main_config = yaml.load(f)

  with open(args.topo_config, 'r') as f:
    topo_config = yaml.load(f)

  return {'main_config': main_config, 'topo_config': topo_config}

def _compute(main_config, topo_config):
  scheme = main_config['scheme']
  solver = None
  logger = logging.getLogger('main')
  if 'output' not in main_config:
    main_config['output'] = 'output.txt'
    logger.warning('No output path provided. Use output.txt as default.')
  if scheme == 'FFC':
    solver = solvers.FFC_Solver(
      main_config=main_config,
      topo_config=topo_config,
      solver_config=None)
  if scheme == 'PCFTF':
    solver = solvers.PCFTF_Solver(
      main_config=main_config,
      topo_config=topo_config,
      solver_config=None)
  if scheme == 'PCFLS':
    solver = solvers.PCFLS_Solver(
      main_config=main_config,
      topo_config=topo_config,
      solver_config=None)
  if scheme == 'PCFLS+':
    solver = solvers.PCFLS_plus_Solver(
      main_config=main_config,
      topo_config=topo_config,
      solver_config=None)
  if solver is None:   
    logger.error('WRONG scheme!')
    return
  mlu, solving_time = solver.compute_mlu(num_link_failure=int(topo_config['attributes']['num_link_failures']))
  logger.info('MLU: %s, solving time: %s seconds' % (mlu, solving_time))
  logger.info('Tunnel reservation is saved in %s' % main_config['output'])

def _main(args, configs):
  main_config = configs['main_config']['main']
  topo_config = configs['topo_config']
  logging.basicConfig(level=main_config['log_level'])
  logger = logging.getLogger('main')

  logger.info("Start.")
  _compute(main_config, topo_config)
  logger.info("Done.")

if __name__ == "__main__":
  args = _parse_args()
  configs = _parse_configs(args)
  _main(args, configs)
