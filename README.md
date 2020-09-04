# PCF

Required library: gurobi, yaml. Please add local gurobi path in pcf/main/run.py if needed.

## How to use:

Execute pcf/main/run.py with a main configuration file and a topology configuration file. For example,

```bash
cd pcf/main
python run.py --main_config ../_config/main.yaml --topo_config ../_config/b4_config.yaml
```

The main configuration file specifies the scheme(FFC, PCFTF, PCFLS, PCFLS+) and the output path. The topology configuration file specifies the topology information including the capacity file, the traffic file and the tunnel file. You will also need to provide the number of link failures you want to design for in the topology configuration file. To use FFC or PCFTF, please also provide the number of tunnels between each pair as a parameter. See pcf/_config/main.yaml and pcf/_config/b4_config.yaml for more detailed format.

## How to generate tunnels:

We provide a tool to generate tunnels and logical sequences in pcf/tools. 

```bash
cd pcf/tools
python generate_tunnels.py --topo_config ../_config/b4_config.yaml --tunnel_type CLS --output_path tunnel.tab
```

Here, you also need to provide a topology configuration file for which you want to generate tunnels. We provide 3 options: 1) physcial: generate 3 physical tunnels per pair; 2) LS: generate physcial tunnels and unconditional logical sequences; 3) CLS: generate physical tunnels and conditional logical sequences. Note that to use this tool, networkx is additionally required. And a local gurobi path also needs to be added if required.

## Configuration file and data file

### main configuration

The main configuration file specifies the scheme(FFC, PCFTF, PCFLS, PCFLS+) and the output path. For example in pcf/_config/main.yaml:

```bash
main:
    log_level: 'INFO'
    output: 'output.txt'     # Output path
#    scheme: 'FFC'
#    scheme: 'PCFTF'
#    scheme: 'PCFLS'
    scheme: 'PCFLS+'
```

### topology configuration

The topology configuration file specifies 

1. Number of physical tunnels between a source-destination pair
2. Number of maximum simultaneous link failures
3. Total number of traffic matrix in the traffic matrix file
4. The index of traffic matrix to be used
5. File paths of capacity file, traffic matrix file, and tunnel file.

For example in pcf/_config/b4_config.yaml:


```bash
name: 'b4'
attributes:
    num_parallel_tunnels: 3   # Number of physical tunnels per source-destination pair
    num_link_failures: 1      # Number of maximum simultaneous link failures 
traffic_matrix:
    num_matrices: 1           # Number of traffic matrix in the traffic file
    tm_index: 0               # Traffic matrix index to be used
data: 
    cap_file: '../_data/b4/b4_capa.tab'               # Capacity file path
    tm_file: '../_data/b4/b4_traffic.tab'             # Traffic matrix file path
    #tunnel_file: '../_data/b4/b4_ffc_tunnel.tab'
    #tunnel_file: '../_data/b4/b4_ls_tunnel.tab'
    tunnel_file: '../_data/b4/b4_cls_tunnel.tab'      # Tunnel file path(including logical sequences)
```


