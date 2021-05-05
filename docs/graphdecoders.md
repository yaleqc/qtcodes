*This file is to take notes while refactoring Graph Decoders. A work in progress...*

# XXZZ Code

## Class Attributes
- `d (int)`: size of surface code
- `T (int)`: number of rounds
- `virtual (dict[str,tuple[float]])`: coordinates of virtual syndrome nodes
- `S (dict[str,nx.Graph])`: dict of syndrome graphs


## Methods
- `correct_readout(readout_string,logical_qubit_value,syndromes)->logical_qubit_value`: corrects logical Z readout string and returns logical_qubit_value as determined by MWPM
-  `corrections(syndromes)->net_flips`: uses the below methods to go from a list of activated syndrome node locations to net_flips on each affected data qubit
   - `make_error_graph(nodes:list[tuple(float)],syndrome_key:str,err_prob:Union[float,None])->error_graph:nx.Graph,paths:`: creates error syndrome subgraph 
     - `_path_degeneracy`: helps make weighted error graph by incorporating deg
   - `matching_graph(error_graph, syndrome_key)->subgraph`: returns matching subgraph
   - `matching(matching_graph,syndrome_key)->filtered_matches`: returns list of relevant matches
   - `calculate_qubit_flips(matches, paths, syndrome_key)->data_qubit_flips_locs`: returns coordinates of data qubits that have undergone flips of syndrome_key type
   -  `net_qubit_flips(flips)->physical_qubit_flips`: returns dict of {data qubit loc:overall flip matrix}

- Graphing
  - `graph_2D(G,edge_lable)`: graphs 2D graphs following certain conventions
  - `graph_2D(G,edge_lable)`: graphs 3D graphs following certain conventions

- Other Helper
  - `_specify_virtual()->None`: populates `virtual`
  - `_make_syndrome_graph()->None`: populates `S`
    - `_populate_syndrome_graph`: helps construct syndrome graphs graphically
    - `_valid_syndrome`: helps construct `S` by checking whether node coordinates are valid
  - `_convert_string_to_nodes(readout_strings)->syndrome_nodes`: convertes readout string to dict of activated syndrome nodes

---

# Repetition Code
- `code`: QEC Code object for which the decoder will be used
- `S` : syndrome nx.Graph describing connectivity between syndrome elements

## Class Attributes

## Methods

- other helper
  - `_separate_string`: string->parsed string
  - `_string2nodes`: string->(nodes, loc of activated sydnrome nodes)
  - `_make_syndrome_graph`:
    - injects all possible Pauli errors in to the circuit for `code`   