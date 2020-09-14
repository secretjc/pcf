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

In general, you can select physical tunnels and logical sequences however you want. We also provide a tool to generate tunnels and logical sequences in pcf/tools as used in PCF paper. 

```bash
cd pcf/tools
python generate_tunnels.py --topo_config ../_config/b4_config.yaml --tunnel_type CLS --output_path tunnel.tab
```

Here, you also need to provide a topology configuration file for which you want to generate tunnels. We provide 3 options: 1) physical: generate 3 physical tunnels per pair; 2) LS: generate physcial tunnels and unconditional logical sequences; 3) CLS: generate physical tunnels and conditional logical sequences. Note that to use this tool, networkx is additionally required. And a local gurobi path also needs to be added if required.

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
    num_matrices: 1           # Number of traffic matrices in the traffic file
    tm_index: 0               # Traffic matrix index to be used
data: 
    cap_file: '../_data/b4/b4_capa.tab'               # Capacity file path
    tm_file: '../_data/b4/b4_traffic.tab'             # Traffic matrix file path
    #tunnel_file: '../_data/b4/b4_ffc_tunnel.tab'
    #tunnel_file: '../_data/b4/b4_ls_tunnel.tab'
    tunnel_file: '../_data/b4/b4_cls_tunnel.tab'      # Tunnel file path(including logical sequences)
```

Note that even though in traffic matrix file, multiple matrices may be provided, only one traffic matrix will be used for one experiment. "tm_index" specifies that particular traffic matrix. See traffic matrix file for more details. 

### capacity file

In the capacity file, each line(in the format of 'i j c') specifies a direct link with capacity of c from i to j. For example in pcf/_data/b4/b4_capa.tab,

```bash
i j capacity
0 1 1
0 2 1
1 0 1
1 4 1
......
```

the first 4 lines specify 4 links 0->1, 0->2, 1->0, 1->4, all with capacity of 1 unit. Note that the unit of capacity is not specified since only the relative values matter. And make sure that the unit of capacity is aligned with the unit of traffic. 

### traffic matrix file

We allow users to provide 1 or more traffic matrices in the traffic matrix file. Different traffic matrices are distinguished by indexes. Each line(in the format of 's t h tm') specifies that in the index h, the demand from s to t is tm. For example in pcf/_data/b4/b4_traffic.tab,

```bash
s t h tm
0 1 0 0.129265573744
0 2 0 0.092814149586
0 3 0 0.063695141693
......
```

the first 3 lines specify the demand from 0->1, 0->2, 0->3 in the traffic matrix 0. Note that, only one traffic matrix will be used for one experiment according to the provided index in the configuration file. 
 
### tunnel file

In the tunnel file, each line either specifies a physical tunnel or a logical sequence. A physical tunnel has 4 entries in a line(source, destination, index, list of edges) and a logical sequence has 5 entries(source, destination, index, list of edges, hints). Note that all physical tunnels need to be specified before any logical sequence. For example in pcf/_data/b4/b4_cls_tunnel.tab,

```bash
s t k edges
7 3 0 7-3
7 3 1 7-6,6-3
7 3 2 7-5,5-2,2-3
......
7 3 0 7-3 no
7 3 1 7-5,5-4,4-3 7-3
......
```

from 7 to 3, there are 3 physical tunnels [7-3], [7-6,6-3] and [7-5,5-2,2-3], and there is one unconditional logical sequence [7-3] and one conditional logical sequence [7-5,5-4,4-3] which is alive when link [7-3] is alive.

Note that you cannot provide any logical sequence if you are using 'FFC' or 'PCFTF' scheme. And if you are using 'PCFLS+' scheme, for every logical sequence entry, two logical sequences will be generated in the model, one is active when the condition is satisfied and the other is active when the condition is not satisfied. 
