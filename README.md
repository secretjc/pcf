# pcf

Required library: gurobi, yaml. Please add local gurobi path in pcf/main/run.py if needed.

## How to use:

Execute pcf/main/run.py with a main configuration file and topology configuration file. For example,

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
