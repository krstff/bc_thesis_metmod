import networkx as nx
import numpy as np
import pandas as pd
import model_polisher as mp
import matplotlib.pyplot as plt
import cobra, libsbml
import json, os, re, escher
from cobra.io import save_json_model, write_sbml_model
from efmtool import calculate_efms

MET_SHEET = "Metabolite List"
RXN_SHEET = "Reaction List"

MET_ID = "ID"
MET_FORMULA = "Formula"
MET_NAME = "Description"
MET_COMPARTMENT = "Compartment"
MET_CHARGE = "Charge"

RXN_ID = "ID"
RXN_NAME = "Description"
RXN_UPPER_BOUND = "Upper bound"
RXN_LOWER_BOUND = "Lower bound"
RXN_OBJECTIVE = "Objective"
RXN_EQUATION = "Reaction"
RXN_SUBSYSTEM = "Subsystem"
RXN_GPR = "GPR"

DEFAULT_COMPARTMENT = "c"
DEFAULT_UPPER_BOUND = 1000
DEFAULT_LOWER_BOUND = -1000

def cobra_to_xlsx(model: cobra.Model, output_path: str) -> None:
    """
    Converts a COBRA model to an Excel file with two sheets: one for reactions and one for metabolites.
    
    Parameters:
    model (cobra.Model): The COBRA model to be converted.
    output_path (str): The path where the Excel file will be saved.

    Returns:
    None
    """
    # Extract reactions
    reaction_data = []
    for rxn in model.reactions:
        reaction_data.append({
            "ID": rxn.id,
            "Reaction": rxn.reaction,
            "GPR": rxn.gene_reaction_rule,
            "Lower bound": rxn.lower_bound,
            "Upper bound": rxn.upper_bound,
            "Objective": rxn.objective_coefficient,
            "Confidence Score": 0,
            "Subsystem": rxn.subsystem,
            "Description": rxn.name
        })

    reactions_df = pd.DataFrame(reaction_data)

    # Extract metabolites
    metabolite_data = []
    for met in model.metabolites:
        metabolite_data.append({
            "ID": met.id,
            "Formula": met.formula,
            "Description": met.name,
            "Compartment": met.compartment,
            "Charge": met.charge
        })

    metabolites_df = pd.DataFrame(metabolite_data)

    # Save to Excel
    with pd.ExcelWriter(output_path) as writer:
        reactions_df.to_excel(writer, sheet_name=RXN_SHEET, index=False)
        metabolites_df.to_excel(writer, sheet_name=MET_SHEET, index=False)

    print(f"Model saved to {output_path}")

def polish_model(model_path: str, polished_path: str, print_diff: bool=False) -> None:
    """Polishes a model using ModelPolisher and saves the polished model to a file.

    Parameters:
    model_path (str): Path to the input SBML model file.
    polished_path (str): Path to save the polished SBML model file.
    print_diff (bool): If True, prints the differences between the original and polished model.

    Returns:
    None
    """
    # Annotation can be further configured - https://github.com/draeger-lab/MPClient/blob/master/examples/default-request-config.json
    polished_result = mp.polish_model_file(model_path)

    polished_model = polished_result["polished_document"]
    success = libsbml.SBMLWriter().writeSBMLToFile(polished_model, polished_path)

    if success:
        print(f"Polished SBML model successfully saved in {polished_path}")
    else:
        print("Failed to save the polished SBML model.")
    if print_diff:
        polished_result["diff"]

def escher_build(model_path: str, map_path: str, online_map: bool=False, highlight_missing: bool=False, reaction_data: dict=None, gene_data: dict=None):
    """Builds an Escher visualisation from a COBRA model and a map file.
    
    Parameters:
    model_path (str): Path to the SBML model file.
    map_path (str): Path to the Escher map file (JSON format).
    online_map (bool): If True, uses an online Escher map.
    highlight_missing (bool): If True, highlights missing reactions in the map.
    reaction_data (dict): Additional data for reactions.
    gene_data (dict): Additional data for genes.
    
    Returns:
    escher.Builder: An Escher map builder object.
    """
    builder = escher.Builder(
        model = cobra.io.read_sbml_model(model_path),
        highlight_missing = highlight_missing,
        reaction_data = reaction_data,
        gene_data = gene_data
    )
    if online_map:
        builder.map_name = map_path
    else:
        builder.map_json = map_path
    return builder


def create_metabolic_graph(model: cobra.Model) -> nx.DiGraph:
    """Creates a directed bipartite graph from a COBRA model.
    Every metabolite and enzyme is a node, their relationships form edges.
    Bipartite graph as there are no connections between enzymes or between metabolites.

    Parameters:
    model (cobra.Model): The COBRA model to be converted into a graph.
    
    Returns:
    nx.DiGraph: A directed bipartite graph representing the model.
    """
    G = nx.DiGraph()

    # Add metabolites as nodes
    for metabolite in model.metabolites:
        G.add_node(
            metabolite.id,
            type="metabolite",
            name=metabolite.name,
            formula=metabolite.formula,
            compartment=metabolite.compartment,
        )

    for reaction in model.reactions:
        # Add enzyme node
        enzyme_node = f"{reaction.id}"
        G.add_node(enzyme_node, type="enzyme", name=reaction.name)
        
        # Add edges from reactants to the reaction
        for reactant in reaction.reactants:
            G.add_edge(reactant.id, enzyme_node, type="reactant", reaction_id=reaction.id)
        
        # Add edges from the reaction to products
        for product in reaction.products:
            G.add_edge(enzyme_node, product.id, type="product", reaction_id=reaction.id)
        
        # Create additional edges for reversible reactions
        if reaction.reversibility:
            for product in reaction.products:
                G.add_edge(product.id, enzyme_node, type="reverse_product", reaction_id=reaction.id)
            for reactant in reaction.reactants:
                G.add_edge(enzyme_node, reactant.id, type="reverse_reactant", reaction_id=reaction.id)

    return G

def xlsx_to_cobra(file_path: str, model_name: str) -> cobra.Model:
    """Converts an Excel file containing metabolic data into a COBRA model.
    The Excel file should contain two sheets: one for reactions and one for metabolites.
    
    Parameters:
    file_path (str): Path to the Excel file.
    model_name (str): Name of the COBRA model to be created.
    
    Returns:
    cobra.Model: A COBRA model object created from the Excel data.
    """
    xls = pd.ExcelFile(file_path)

    # Load pages
    reactions_df = pd.read_excel(xls, sheet_name=RXN_SHEET)
    metabolites_df = pd.read_excel(xls, sheet_name=MET_SHEET)

    # Create empty COBRA model
    model = cobra.Model(model_name)

    # Add metabolites to the model
    for _, row in metabolites_df.iterrows():
        met = cobra.Metabolite(
            id=row[MET_ID],
            formula=row[MET_FORMULA],
            name=row[MET_NAME],
            compartment=row[MET_COMPARTMENT] if not pd.isna(row[MET_COMPARTMENT]) else DEFAULT_COMPARTMENT,
            charge=row[MET_CHARGE],
        )
        model.add_metabolites([met])

    # Add reactions to the model
    for _, row in reactions_df.iterrows():
        reaction = cobra.Reaction(id=row[RXN_ID])
        reaction.name = row[RXN_NAME]
        reaction.lower_bound = row[RXN_LOWER_BOUND] if not pd.isna(row[RXN_LOWER_BOUND]) else DEFAULT_LOWER_BOUND
        reaction.upper_bound = row[RXN_UPPER_BOUND] if not pd.isna(row[RXN_UPPER_BOUND]) else DEFAULT_UPPER_BOUND
        if not pd.isna(row[RXN_GPR]):
            reaction.gene_reaction_rule=row[RXN_GPR]
        reaction.subsytem=row[RXN_SUBSYSTEM]

        model.add_reactions([reaction])
        model.reactions.get_by_id(row[RXN_ID]).build_reaction_from_string(row[RXN_EQUATION]) # This function has to have access to model - therefore the reaction needs to be already present in model
        model.reactions.get_by_id(row[RXN_ID]).objective_coefficient=row[RXN_OBJECTIVE] if not pd.isna(row[RXN_OBJECTIVE]) else 0 # Same here

    return(model)

def print_graph_nodes_names(G: nx.Graph) -> None:
    print(nx.get_node_attributes(G, "name"))

def find_graph_names_with(G: nx.DiGraph, str: str) -> list:
    """Returns a list of names of nodes that contain the given string
    in their name"""
    return [name for name in nx.get_node_attributes(G, "name") if str in name]

def print_neighbors(G: nx.Graph, node: str):
    for n in G.neighbors(node):
        print(n)

def get_id_by_name(G: nx.Graph, name: str) -> str:
    """Returns the id of a node with the given name.
    If there are multiple nodes with the same name, returns the first one.
    """
    names = nx.get_node_attributes(G, "name")
    for k, v in names.items():
        if v == name:
            return k
    raise Exception(f"Name {name} not found")

def get_name_by_id(G: nx.Graph, id: str) -> str:
    """Returns the name of a node with the given id."""
    names = nx.get_node_attributes(G, "name")
    return names[id]

def print_degree_distribution(G: nx.Graph) -> None:
    degrees = [degree for _, degree in G.degree()]
    plt.hist(degrees, bins=range(min(degrees), max(degrees) + 1), edgecolor="black")
    plt.title("Node Degree Distribution")
    plt.xlabel("Degree")
    plt.ylabel("Frequency")
    plt.show()

def print_degree_sorted(G: nx.Graph) -> None:
    node_degrees = [(node, G.in_degree(node) + G.out_degree(node)) for node in G.nodes()]
    sorted_nodes_reverse = sorted(node_degrees, key=lambda x: x[1], reverse=True)
    sorted_nodes_non_reverse = sorted(node_degrees, key=lambda x: x[1])
    print("Highest degree nodes:")
    print(sorted_nodes_reverse)
    print("Lowest degree nodes:")
    print(sorted_nodes_non_reverse)

def print_strongly_connected(G: nx.Graph) -> None:
    """Prints all strongly connected components containing more than one node."""
    nx.number_strongly_connected_components(G)
    sorted_c = sorted(list(nx.strongly_connected_components(G)), key=len)
    for comp in sorted_c:
        # Ignore single node components
        if len(comp) > 1:
            print(len(comp), comp)

def find_paths(G: nx.Graph, source: str, target: str, nodes: list, count=5, print_path=False) -> tuple:
    """Finds paths in a directed graph from source to target, including all nodes.

    Parameters:
    G (nx.Graph): The graph.
    source (str): The starting node.
    target (str): The target node.
    nodes (list): List of nodes that must be included in the path.
    count (int): Number of paths to find.
    print_path (bool): If True, prints the paths.

    Returns:
    tuple: Tuple of paths found. First one uses IDs, second one uses names.
    """
    paths = []
    named_paths = []

    for path in nx.shortest_simple_paths(G, source, target):
        skip = False
        for node in nodes:
            if node not in path:
                skip = True
                break

        if skip:
            continue

        count -= 1
        paths.append(path)
        named_path = [get_name_by_id(G, node_id) for node_id in path]
        named_paths.append(named_path)
        
        if print_path:
            print(named_path)

        if count == 0:
            break
    
    return paths, named_paths

def delete_nodes_by_id(G: nx.Graph, ids: list) -> nx.Graph:
    for id_ in ids:
        G.remove_node(id_)
    return G

def prune_graph(G: nx.Graph, threshold: int, to_keep: list) -> nx.Graph:
    """Prunes a graph by removing nodes with degree higher than a threshold.
    Nodes in the to_keep list will not be removed, even if they exceed the threshold.

    Parameters:
    G (nx.Graph): The graph to be pruned.
    threshold (int): The degree threshold for pruning.
    to_keep (list): List of nodes that should not be removed.
    
    Returns:
    nx.Graph: The pruned graph.
    """
    G_pruned = G.copy()

    high_degree_nodes = [node for node, degree in G.degree() if (degree > threshold and node not in to_keep) or "UDP" in G.nodes[node].get("name", "")]
    to_delete_named = [get_name_by_id(G, node) for node in high_degree_nodes]

    G_pruned.remove_nodes_from(high_degree_nodes)

    print(f"Removed {len(high_degree_nodes)} nodes with degree greater than {threshold}.")
    print(f"Removed the following nodes: {to_delete_named}")
    print(f"Remaining nodes: {len(G.nodes())}")
    return G_pruned

def sync_xlsx(input_file: str, output_name: str, save_json: bool=True) -> None:
    """Coverts xlsx file to json and SBML.

    Parameters:
    input_file (str): Path to the input xlsx file.
    output_name (str): Name for the output files (without extension).
    save_json (bool): If True saves the model to JSON.

    Returns:
    None
    """
    # coverts xlsx file to json and SBML
    model = xlsx_to_cobra(input_file, output_name)
    write_sbml_model(model,  output_name + ".xml")
    if save_json:
        save_json_model(model, output_name + ".json")
    print("Saved")

def sbml_to_json(model_path: str, output_path: str) -> None:
    model = cobra.io.read_sbml_model(model_path)
    save_json_model(model, output_path)

def json_to_sbml(model_path: str, output_path: str) -> None:
    model = cobra.io.json.load_json_model(model_path)
    cobra.io.write_sbml_model(model, output_path)

def mini_model_prune(escher_map, model_path, output_name):
    """Prunes a COBRA model based on an Escher map.
    Deletes all reactions and metabolites not present in the Escher map and saves the pruned model in XML, JSON, and XLSX formats.
    
    Parameters:
    escher_map (str): Path to the Escher map file (JSON format).
    model_path (str): Path to the SBML model file.
    output_name (str): Name for the output files (without extension).
    
    Returns:
    None
    """
    model = cobra.io.read_sbml_model(model_path)

    with open(escher_map, "r") as file:
        json_data = json.load(file)

    map_reactions = set(reaction_info["bigg_id"] for reaction_info in json_data[1].get("reactions", {}).values())

    reactions_to_remove = [rxn.id for rxn in model.reactions if rxn.id not in map_reactions and "Biomass" not in rxn.name]
    model.remove_reactions(reactions_to_remove)

    orphan_metabolites = [met for met in model.metabolites if len(met.reactions) == 0]
    model.remove_metabolites(orphan_metabolites)

    # Ensure models directory exists
    os.makedirs("models", exist_ok=True)
    base = os.path.join("models", output_name)
    # Save in xml, json and xlsx
    cobra.io.write_sbml_model(model, base + ".xml")
    cobra.io.save_json_model(model, base + ".json")
    cobra_to_xlsx(model, base + ".xlsx")

    print(f"Filtered model saved as {base}")
    print(f"Removed {len(reactions_to_remove)} reactions.")
    print(f"Removed {len(orphan_metabolites)} metabolites.")

def add_objective_function(formula: str, model: cobra.Model, name: str="Objective_function", id: str="obj_rxn") -> None:
    """Constructs a reaction from formula string, adds it to a model and sets it as objective."""
    reaction = cobra.Reaction(id)
    reaction.name = name
    reaction.lower_bound = 0
    reaction.upper_bound = 1000
    model.add_reactions([reaction])

    reaction.build_reaction_from_string(formula)
    model.objective = reaction

def decode_sbml(str: str) -> str:
    return str.replace("__45__", "-").replace("__46__", ".").replace("__95__", "_")

def encode_sbml(str: str) -> str:
    return str.replace("_", "__95__").replace("+", "__43__").replace(",", "__44__").replace("-", "__45__").replace(".", "__46__")

def extract_bigg_id(annotation, keyword):
    """Helper function.
    Extracts BiGG ID from the annotation of a reaction or metabolite.
    """
    # note that ChatGPT was used to discuss this function heavily
    if annotation is None:
        return None

    for i in range(annotation.getNumChildren()):
        child = annotation.getChild(i)
        if child.getName() == "RDF":
            for j in range(child.getNumChildren()):
                desc = child.getChild(j)
                if desc.getName() == "Description":
                    for k in range(desc.getNumChildren()):
                        bag_wrapper = desc.getChild(k)
                        if bag_wrapper.getName() == "is":
                            bag = bag_wrapper.getChild(0)
                            for l in range(bag.getNumChildren()):
                                li = bag.getChild(l)
                                attrs = li.getAttributes()
                                for idx in range(attrs.getLength()):
                                    name = attrs.getName(idx)
                                    value = attrs.getValue(idx)
                                    uri = attrs.getURI(idx)
                                    if name == "resource" and "rdf" in uri and keyword in value:
                                        return value.split("/")[-1]
    return None

def convert_gpr_to_string(assoc: libsbml.GeneProductAssociation) -> str:
    """Helper function.
    Converts a gene product association to a string.
    """
    if assoc is None:
        return ""
    if assoc.getTypeCode() == libsbml.AND_ASSOCIATION:
        return "(" + " and ".join([convert_gpr_to_string(assoc.getAssociation(i)) for i in range(assoc.getNumAssociations())]) + ")"
    elif assoc.getTypeCode() == libsbml.OR_ASSOCIATION:
        return "(" + " or ".join([convert_gpr_to_string(assoc.getAssociation(i)) for i in range(assoc.getNumAssociations())]) + ")"
    elif isinstance(assoc, libsbml.GeneProductRef):
        return assoc.getGeneProduct()
    return ""

def build_bigg_id_maps(sbml_path: str) -> tuple:
    """Helper function.
    Builds a mapping of BiGG IDs to SBML IDs for reactions and metabolites.
    """
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(sbml_path)
    model = doc.getModel()

    reaction_id_map = {}
    metabolite_id_map = {}

    for reaction in model.getListOfReactions():
        bigg_id = extract_bigg_id(reaction.getAnnotation(), "bigg.reaction")
        if bigg_id:
            reaction_id_map[bigg_id] = reaction.getId()

    for species in model.getListOfSpecies():
        bigg_id = extract_bigg_id(species.getAnnotation(), "bigg.metabolite")
        if bigg_id:
            metabolite_id_map[bigg_id] = species.getId()

    return reaction_id_map, metabolite_id_map

def escher_map_sync(escher_map_path: str, model_path: str) -> None:
    """Synchronizes an Escher map with a COBRA model by replacing metabolite and reaction names in the map with those from the model.
    The function reads the Escher map and the COBRA model, builds a mapping of BiGG IDs to SBML IDs, and updates the Escher map accordingly.
    The final map is saved in maps/ directory.
    
    Parameters:
    escher_map_path (str): Path to the Escher map file (JSON format).
    model_path (str): Path to the SBML model file.
    
    Returns:
    None
    """
    def replace_markers_in_file(filename, output_filename):
        with open(filename, "r", encoding="utf-8") as f:
            content = f.read()

        content = decode_sbml(content)

        with open(output_filename, "w", encoding="utf-8") as f:
            f.write(content)
    
    def replace_values(obj, replace_dict):
        if isinstance(obj, dict):
            return {k: replace_values(v, replace_dict) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [replace_values(elem, replace_dict) for elem in obj]
        elif isinstance(obj, str):
            for key, value in replace_dict.items():
                # It seems that Escher automatically removes these prefixes, therefore
                # we can't add them to the map
                if "R_" in value or "M_" in value or "G_" in value:
                    value = value[2:]
                obj = obj.replace(key, value)
            return obj
        else:
            return obj
    
    # SBML changes some characters (like "-") to a standardized naming scheme,
    # sometimes this can produce errors https://github.com/opencobra/cobrapy/issues/936
    replace_markers_in_file(model_path, "temp.xml")

    # a mapping of reaction and metabolite names is created
    rxn_map, met_map = build_bigg_id_maps("temp.xml")

    with open(escher_map_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    updated_data = replace_values(data, (met_map | rxn_map))

    # Ensure maps directory exists
    os.makedirs("maps", exist_ok=True)
    output_map_path = os.path.join("maps", f"renamed_{os.path.basename(escher_map_path)}")
    with open(output_map_path, "w", encoding="utf-8") as f:
        json.dump(updated_data, f, indent=2)
    
    os.remove("temp.xml")
    print(f"Updated map saved as {output_map_path}")

def build_gene_map(sbml_path: str) -> dict:
    """Helper function.
    Builds a mapping of reactions to gene product associations.
    """
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(sbml_path)
    model = doc.getModel()

    rxn_gpr_map = {}

    for reaction in model.getListOfReactions():
        fbc_rxn_plugin = reaction.getPlugin("fbc")
        assoc = fbc_rxn_plugin.getGeneProductAssociation()
        if assoc:
            gpr_str = convert_gpr_to_string(assoc.getAssociation())
            rxn_gpr_map[decode_sbml(reaction.getId().removeprefix("R_"))] = gpr_str.removeprefix("G_")

    return rxn_gpr_map

def replace_genes_in_escher_map(escher_map_path: str, gene_mapping: dict, output_path: str, to_decode_sbml: bool) -> None:
    """Helper function.
    Replaces gene reaction rules in an Escher map with a mapping from a COBRA model.
    """
    with open(escher_map_path, "r", encoding="utf-8") as f:
        escher_raw = json.load(f)

    escher_data = escher_raw[1]

    for _, rxn_data in escher_data.get("reactions", {}).items():
        rxn_id = rxn_data["bigg_id"]
        if rxn_id in gene_mapping:
            genes = gene_mapping[rxn_id]
            if to_decode_sbml:
                genes = decode_sbml(genes)

            rxn_data["gene_reaction_rule"] = genes
            gene_ids = set(re.findall(r"\b\w+\b", genes))
            logic_keywords = {"and", "or"}
            gene_ids = [g for g in gene_ids if g.lower() not in logic_keywords]

            # Build gene list
            rxn_data["genes"] = [
                {
                    "name": f"G_{gene_id}",
                    "bigg_id": gene_id
                }
                for gene_id in sorted(gene_ids)
            ]

    escher_raw[1] = escher_data
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(escher_raw, f, indent=2)

def sync_map_genes(model_path: str, escher_map_path: str, output_path: str="updated_map.json", to_decode_sbml: bool=True) -> None:
    """Synchronizes an Escher map with a COBRA model by replacing gene reaction rules in the map with those from the model.
    The function reads the Escher map and the COBRA model, builds a mapping of gene reaction rules, and updates the Escher map accordingly.
    The final map is saved in maps/ directory.
    
    Parameters:
    model_path (str): Path to the SBML model file.
    escher_map_path (str): Path to the Escher map file (JSON format).
    output_path (str): Path to save the updated Escher map.
    decode_sbml (bool): If True, decodes SBML IDs in the mapping.
    
    Returns:
    None
    """
    gene_map = build_gene_map(model_path)
    rxn_map, _ = build_bigg_id_maps(model_path)
    
    id_to_gene = {}
    for bigg_id, model_id in rxn_map.items():
        model_id = decode_sbml(model_id)
        gene = gene_map.get(model_id)
        if gene:
            id_to_gene[bigg_id] = gene
    # Ensure maps directory exists
    os.makedirs("maps", exist_ok=True)
    output_path_final = output_path
    if not os.path.dirname(output_path):
        output_path_final = os.path.join("maps", output_path)
    replace_genes_in_escher_map(escher_map_path, gene_map | id_to_gene, output_path_final, to_decode_sbml)
    print(f"Gene rules updated and saved to {output_path_final}")

def filter_reaction(reaction: str, model: cobra.Model, arrow_symbol: str="→") -> str:
    """Deletes metabolites from a reaction string that are not in the model."""
    reactants, products = reaction.split(arrow_symbol)
    
    def filter_metabolites(segment: str) -> str:
        metabolite_pattern = re.findall(r"([\d\.]+)\s+(\w+)", segment)  # Extract (coefficient, metabolite) pairs
        filtered_metabolites = [f"{coef} {met}" for coef, met in metabolite_pattern if met in model.metabolites]
        return " + ".join(filtered_metabolites)
    
    filtered_reactants = filter_metabolites(reactants)
    filtered_products = filter_metabolites(products)

    return f"{filtered_reactants} --> {filtered_products}" if filtered_reactants or filtered_products else ""

def filter_reaction_bigg(reaction: str, model_path: str, arrow_symbol: str="→"):
    """Deletes metabolites from a reaction string that are not in the model.
    The function is similar to filter_reaction, but it uses BiGG IDs instead of SBML IDs and renames all metabolites to use SBML IDs.
    """
    _, met_map = build_bigg_id_maps(model_path)
    model = cobra.io.read_sbml_model(model_path)

    reactants, products = reaction.split(arrow_symbol)

    def remove_compartment_suffix(met_id: str) -> str:
        return met_id.rsplit("_", 1)[0]

    def filter_metabolites(segment: str) -> str:
        metabolite_pattern = re.findall(r"([\d\.]+)\s+(\w+)", segment)
        filtered_metabolites = []
        for coef, met in metabolite_pattern:
            met_core = remove_compartment_suffix(met)
            if met_core in met_map and met_map[met_core] in model.metabolites:
                filtered_metabolites.append(f"{coef} {met_map[met_core]}")
        return " + ".join(filtered_metabolites)

    filtered_reactants = filter_metabolites(reactants)
    filtered_products = filter_metabolites(products)

    return f"{filtered_reactants} --> {filtered_products}" if filtered_reactants or filtered_products else ""

def calculate_efm_matrix(model: cobra.Model, log: bool=False):
    """Calculates the elementary flux modes (EFMs) for a given COBRA model using the efmtool library.
    
    Parameters:
    model (cobra.Model): The COBRA model for which to calculate EFMs.
    log (bool): If True, logs the output to the console.
    
    Returns:
    np.ndarray: A matrix where each column represents an EFM and each row represents a reaction.
    """
    S = cobra.util.array.create_stoichiometric_matrix(model)

    reversibilities = [1 if rxn.reversibility else 0 for rxn in model.reactions]

    # Extract reaction and metabolite names
    reaction_names = [rxn.id for rxn in model.reactions]
    metabolite_names = [met.id for met in model.metabolites]

    level = "WARNING"
    if log:
        level = "INFO"
    
    efms = calculate_efms(
        stoichiometry=S,
        reversibilities=reversibilities,
        reaction_names=reaction_names,
        metabolite_names=metabolite_names,
        options={"kind": "stoichiometry",
                "arithmetic": "double",
                "zero": "1e-10",
                "compression": "default",
                "log": "console",
                "level": level, # Changed from default in order to suppress logger
                "maxthreads": "-1",
                "normalize": "min",
                "adjacency-method": "pattern-tree-minzero",
                "rowordering": "MostZerosOrAbsLexMin"}
    )
    print(f"Matrix size: {efms.shape}")
    return efms

def find_smallest_efm(efm_matrix: np.ndarray) -> int:
    """Finds the index of the smallest elementary flux mode (EFM) in the matrix.
    Use efm = efm_matrix[:, idx]
    """
    num_zeros = np.sum(efm_matrix == 0, axis=0)
    idx = np.argmax(num_zeros)
    num_nonzeros = efm_matrix.shape[0] - num_zeros[idx]
    print(f"EFM at index {idx} has {num_nonzeros} non-zero (active) reactions.")
    return idx


def efm_for_escher(model: cobra.Model, efm_matrix: np.ndarray, save_to_file: bool=False, output_path: str="efm.json") -> dict:
    """Saves a single elementary flux mode (EFM) as a JSON file for Escher visualisation.

    Parameters:
    model (cobra.Model): The COBRA model.
    efm_matrix (np.ndarray): The EFM matrix.
    save_to_file (bool): If True, saves the EFM to a file.
    output_path (str): The path to save the EFM file.

    Returns:
    dict: A dictionary mapping reaction names to their corresponding flux values in the EFM.
    """
    reaction_names = [rxn.id for rxn in model.reactions]

    reaction_to_flux = dict(zip(reaction_names, efm_matrix))
    if save_to_file:
        with open(output_path, "w") as json_file:
            json.dump(reaction_to_flux, json_file, separators=(",", ":"), 
                    default=lambda x: f"{x:.6e}")
    return reaction_to_flux

def filter_efms(efm_matrix: np.ndarray, rxn_list: list, model: cobra.Model) -> np.ndarray:
    """Filters the EFM matrix to include only reactions listed in rxn_list.

    Parameters:
    efm_matrix (np.ndarray): The EFM matrix.
    rxn_list (list): A list of reaction IDs to keep in the EFM matrix.
    model (cobra.Model): The COBRA model.

    Returns:
    np.ndarray: The filtered EFM matrix.
    """
    reaction_names = [rxn.id for rxn in model.reactions]
    reaction_names_map = {rxn: idx for idx, rxn in enumerate(reaction_names)}
    matching_indices = []

    for i in range(efm_matrix.shape[1]):
        column = efm_matrix[:, i]

        if all(column[reaction_names_map[rxn]] != 0 for rxn in rxn_list):
            matching_indices.append(i)

    return efm_matrix[:, matching_indices]

def combine_all_efms(efm_matrix: np.ndarray) -> np.ndarray:
    """Combines all elementary flux modes (EFMs) into a single vector."""
    return np.sum(efm_matrix, axis=1)

def find_reaction_not_in_efms(efm_matrix: np.ndarray, model: cobra.Model) -> list:
    # Finds reactions that are not present in any elementary flux mode (EFM).
    reaction_names = [rxn.id for rxn in model.reactions]

    vector = combine_all_efms(efm_matrix)
    zero_indices = np.where(vector == 0)[0]
    zero_reactions = [reaction_names[i] for i in zero_indices]
    return zero_reactions

def filter_gene_csv(input_file: str, col1_idx: int, col2_idx: int, output_file: str="filtered_file.csv",
                     delim: str=",", prefix: str="", num_prefix: str="", encode_sbml: bool=False):
    """Filters a CSV file to keep only two specified columns and renames them.
    The first column is used as an ID, and the other two columns are renamed to "re_x" and "re_y".
    https://escher.readthedocs.io/en/latest/getting_started.html#creating-data-files-as-csv-and-json
    The function also encodes the IDs based on the provided prefix and number prefix.
    If GPR in model is in form XX_12345_YY12345 use prefix for XX and num_prefix for YY
    
    Parameters:
    input_file (str): Path to the input CSV file.
    col1_idx (int): Index of the first column to keep.
    col2_idx (int): Index of the second column to keep.
    output_file (str): Path to save the filtered CSV file.
    delim (str): Delimiter used in the input CSV file.
    prefix (str): Prefix to add to the IDs.
    num_prefix (str): Number prefix to add to the IDs.
    encode_sbml (bool): If True, encodes the IDs using SBML encoding.
    
    Returns:
    list: A list of two dictionaries, where the first dictionary maps IDs to "re_x" values,
          and the second dictionary maps IDs to "re_y" values.
    """
    df = pd.read_csv(input_file, delimiter=delim)
    
    id_column = df.columns[0]
    col1 = df.columns[col1_idx]
    col2 = df.columns[col2_idx]
    
    # Create a new DataFrame with the selected columns

    new_df = df[[id_column, col1, col2]].copy()
    new_df.columns = ["id", "re_x", "re_y"]

    # The model saves info about gene data in "geneProductAssociation" attribute
    # Therefore the same name has to be in the gene data csv
    if encode_sbml:
        new_df["id"] = new_df["id"].apply(lambda x: encode_sbml(x))
        new_df["id"] = new_df["id"].apply(lambda x: f"{prefix}{x}" if num_prefix in x else f"{prefix}{x.replace('__95__', encode_sbml(num_prefix))}")
    else:
        new_df["id"] = new_df["id"].apply(lambda x: f"{prefix}{x}" if num_prefix in x else f"{prefix}{x.replace('_', num_prefix)}")


    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    new_df.to_csv(output_file, index=False)
    print(f"Filtered CSV saved as {output_file}")

    dict1 = dict(zip(new_df["id"], new_df["re_x"]))
    dict2 = dict(zip(new_df["id"], new_df["re_y"]))

    result = [dict1, dict2]

    return result