import glob, json, os, subprocess, tempfile, argparse, math
import networkx as nx
import numpy as np
import pandas as pd
from pyvis.network import Network


# Return species from identifier tag for each PDB entry.
def get_species(protein_id, species_id):
    for species in species_id:
        if protein_id.startswith(
            species
        ):  # e.g. "test_sponsalis_29" -> "Conus sponsalis structure no. 29"
            return species
    return "Unknown"


# Load Foldseek results from file, looking for query, targets and E-values.
def load_foldseek(files, evalue_threshold=1.0e-3, species_id=[], chunk_size=100000):
    chunks = []  # List to accumulate chunks of data
    dtypes = {0: str, 1: str, 2: float}
    col_names = [
        "query",
        "target",
        "evalue",
    ]  # Foldseek --output-format should contain these columns

    for file in files:
        for chunk in pd.read_csv(
            file,
            sep="\t",
            header=None,
            names=col_names,
            chunksize=chunk_size,
            dtype=dtypes,
        ):
            # Filter and discard unwanted data
            chunk = chunk[
                chunk["evalue"] <= evalue_threshold
            ]  # Filter by E-value threshold (set to 1.0e-3)
            chunk = chunk[chunk["query"] != chunk["target"]]
            # Process species for the chunk.
            chunk["query_species"] = chunk["query"].apply(
                lambda x: get_species(x, species_id)
            )
            chunk["target_species"] = chunk["target"].apply(
                lambda x: get_species(x, species_id)
            )
            chunks.append(chunk)
    if chunks:
        df = pd.concat(
            chunks, ignore_index=True
        )  # Concatenate all chunks into a single DataFrame
    else:
        df = pd.DataFrame(columns=col_names)
    return df


def build_graph(df):
    G = nx.Graph()  # Create an empty graph, G
    for _, row in df.iterrows():
        query, target, evalue = row["query"], row["target"], row["evalue"]
        if query not in G:  # Ensure nodes are created only once
            G.add_node(query, species=row["query_species"])
        if target not in G:
            G.add_node(target, species=row["target_species"])
        G.add_edge(query, target, weight=evalue)
    return G


def superpose(query_pdb, target_pdb):
    # Run PyMOL to superpose two PDB structures and return the aligned PDB content.
    try:
        with tempfile.NamedTemporaryFile(
            delete=False, suffix=".pdb"
        ) as fq:  # Create a temporary file for the query PDB
            fq.write(query_pdb.encode("utf-8"))
            query_file = fq.name
        with tempfile.NamedTemporaryFile(
            delete=False, suffix=".pdb"
        ) as ft:  # Create a temporary file for the target PDB
            ft.write(target_pdb.encode("utf-8"))
            target_file = ft.name
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as fo:
            out_file = fo.name

        # PyMOL script to perform the superposition
        pymol_script = f"""
            load {query_file}, query
            load {target_file}, target
            super query, target
            save {out_file}, query
            quit
        """
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pml") as fs:
            fs.write(pymol_script.encode("utf-8"))
            script_file = fs.name

        # Run PyMOL with the script
        subprocess.run(
            ["pymol", "-qc", script_file],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
        )
        with open(out_file, "r") as f_out:
            result = f_out.read()
    except Exception as e:
        print("Error running PyMOL superposition:", e)
        result = query_pdb
    finally:
        for f in [query_file, target_file, out_file, script_file]:
            try:
                os.remove(f)
            except Exception:
                pass
    return result


# Load classification file and return mapping
def load_AFDBhits(file_path):
    df_af = pd.read_csv(
        file_path, sep="\t", header=None, names=["query", "target", "evalue"]
    )
    afdb_classification = {}
    for _, row in df_af.iterrows():
        query = row["query"]  # e.g. "test_sponsalis_29"
        entry = {
            "target": row["target"],
            "evalue": row["evalue"],
        }  # e.g. "AF-{acc.id}-F1-model_v4"
        afdb_classification.setdefault(query, []).append(entry)
    return afdb_classification


# Define mappings for VenomZone and Conotoxin classifications
toxinfold_map = {
    "Kunitz.tsv": "Venom Kunitz-type",
    "Flavin_monoamine_oxidase.tsv": "Flavin monoamine oxidase",
    "Endothelin-sarafotoxin.tsv": "Endothelin/sarafotoxin",
    "Arthropod_hormone.tsv": "Arhtropod hormone",
    "Arthropod_dermonecrotic_toxin.tsv": "Arhtropod dermonecrotic toxin",
    "AVIT_prokineticin.tsv": "AVIT (prokineticin)",
    "Bradykinin-potentiating_peptide.tsv": "Bradykinin-potentiating peptide (BPP)",
    "Bradykinin-related_peptide.tsv": "Bradykinin-related peptide",
    "Cathelicidin.tsv": "Cathelicidin",
    "Complement_C3.tsv": "Complement C3 homolog",
    "CRISP.tsv": "CRISP",
    "Crotamine-myotoxin.tsv": "Crotamine-myotoxin",
    "Coninsulin.tsv": "Coninsulin",
    "Cystatin.tsv": "Cystatin",
    "Disintegrin.tsv": "Disintegrin",
    "5_nucleotidase.tsv": "5'-nucleotidase",
    "AB_hydrolase.tsv": "AB hydrolase superfamily",
    "EGF_domain_peptide.tsv": "EGF domain peptide",
    "ICK.tsv": "ICK (Knottin)",
    "Glycosyl_hydrolase_56.tsv": "Glycosyl hyrolase 56 (hyaluronidase)",
    "MetalloproteinaseM12B.tsv": "Metalloproteinase (M12B)",
    "NDBP.tsv": "Scorpion non-disulfide-bridged peptide (NDBP)",
    "NGF_beta.tsv": "NGF-beta",
    "Ohanin-vespryn.tsv": "Ohanin/Vespryn",
    "Nucleotide_pyrophosphatase.tsv": "Nucleotide pyrophosphatase",
    "Natriuretic_peptide.tsv": "Natriuretic peptide",
    "NaChannelAlpha.tsv": "Scorpion Na+ channel toxin, alpha",
    "NaChannelBeta.tsv": "Scorpion Na+ channel toxin, beta",
    "Pacifastin_protease.tsv": "Pacifastin protease inhibitor",
    "True_venom_lectin.tsv": "True venom lectin",
    "ThreeFingerToxin.tsv": "Three-finger toxin",
    "TCTP.tsv": "Translationally controlled tumor protein",
    "pHpG.tsv": "pHpG (metalloproteinase inhibitor)",
    "Scorpktx.tsv": "Scorpion potassium channel toxin",
    "Snaclec.tsv": "Snaclec",
    "Peptidase-S1.tsv": "Peptidase S1 (serine protease)",
    "Phospholipase-A2.tsv": "Phospholipase A2",
    "scrip.tsv": "Cnidaria small cysteine-rich protein",
    "PDGF-VEGF.tsv": "PDGF/VEGF growth factor",
    "Multicopper_oxidase.tsv": "Multicopper oxidase",
    # "Classification.tsv": "Classification" - for further entries
}

superfam_map = {
    "supFamA.tsv": "Superfamily A",
    "supFamB.tsv": "Superfamily B",
    "supFamB2.tsv": "Superfamily B2",
    "supFamB4.tsv": "Superfamily B4",
    "supFamC.tsv": "Superfamily C",
    "supFamD.tsv": "Superfamily D",
    "supFamE.tsv": "Superfamily E",
    "supFamF.tsv": "Superfamily F",
    "supFamH.tsv": "Superfamily H",
    "supFamJ.tsv": "Superfamily J",
    "supFamK.tsv": "Superfamily K",
    "supFamL.tsv": "Superfamily L",
    "supFamL11.tsv": "Superfamily L11",
    "supFamI1.tsv": "Superfamily I1",
    "supFamI2.tsv": "Superfamily I2",
    "supFamI3.tsv": "Superfamily I3",
    "supFamI4.tsv": "Superfamily I4",
    "supFamConikotikot.tsv": "Superfamily Conikotikot",
    "supFamConopressin.tsv": "Superfamily Conopressin",
    "supFamConodipine.tsv": "Superfamily Conodipine",
    "supFamInsulin.tsv": "Superfamily Coninsulin",
    "supFamConorfamide.tsv": "Superfamily Conorfamide",
    "supFamConocap.tsv": "Superfamily Conocap",
    "supFamConkunitzin.tsv": "Superfamily Conkunitzin",
    "supFamConohyal.tsv": "Superfamily Conohyal",
    "supFamConoporin.tsv": "Superfamily Conoporin",
    "supFamM.tsv": "Superfamily M",
    "supFamS.tsv": "Superfamily S",
    "supFamMARFL.tsv": "Superfamily MARFL",
    "supFamMASEG.tsv": "Superfamily MASEG",
    "supFamMEFFR.tsv": "Superfamily MEFFR",
    "supFamMGATL.tsv": "Superfamily MGATL",
    "supFamMGGRF.tsv": "Superfamily MGGRF",
    "supFamMKAVA.tsv": "Superfamily MKAVA",
    "supFamMKFLL.tsv": "Superfamily MKFLL",
    "supFamMKISL.tsv": "Superfamily MKISL",
    "supFamMKIYL.tsv": "Superfamily MKIYL",
    "supFamMLLLL.tsv": "Superfamily MLLLL",
    "supFamMLSML.tsv": "Superfamily MLSML",
    "supFamMMLFM.tsv": "Superfamily MMLFM",
    "supFamMo3964.tsv": "Superfamily Mo3964",
    "supFamMr30.tsv": "Superfamily Mr30",
    "supFamMRFYM.tsv": "Superfamily MRFYM",
    "supFamMSCLY.tsv": "Superfamily MSCLY",
    "supFamMSRLF.tsv": "Superfamily MSRLF",
    "supFamMSRSG.tsv": "Superfamily MSRSG",
    "supFamMWIRK.tsv": "Superfamily MWIRK",
    "supFamN.tsv": "Superfamily N",
    "supFamO1.tsv": "Superfamily O1",
    "supFamO2.tsv": "Superfamily O2",
    "supFamO3.tsv": "Superfamily O3",
    "supFamP.tsv": "Superfamily P",
    "supFamProhormone.tsv": "Superfamily Prohormone 4",
    "supFamQ.tsv": "Superfamily Q",
    "supFamSF-04.tsv": "Superfamily SF-04",
    "supFamSF-mi1.tsv": "Superfamily SF-mi1",
    "supFamSF-mi2.tsv": "Superfamily SF-mi2",
    "supFamT.tsv": "Superfamily T",
    "supFamU.tsv": "Superfamily U",
    "supFamunk1.tsv": "Superfamily unk1",
    "supFamunk2.tsv": "Superfamily unk2",
    "supFamV.tsv": "Superfamily V",
    "supFamY.tsv": "Superfamily Y",
    # "Classification.tsv": "Classification" - for further entries
}


def vzone_map(mapping_dir, file_mapping):
    # Build classification for VenomZone folds
    mapping = {}
    for filename, classification in file_mapping.items():
        filepath = os.path.join(mapping_dir, filename)
        if os.path.exists(filepath):
            with open(filepath, "r") as f:
                for line in f:
                    raw_id = line.strip()
                    if raw_id:
                        # Convert the raw ID to canonical annotation
                        processed_id = f"AF-{raw_id}-F1-model_v4"
                        mapping[processed_id] = classification
        else:
            print(f"Warning: Mapping file {filepath} not found.")
    return mapping


def cono_map(mapping_dir, file_mapping):
    # Build classification for Conotoxin folds
    mapping_cono = {}
    for filename, classification in file_mapping.items():
        filepath = os.path.join(mapping_dir, filename)
        if not os.path.exists(filepath):
            print(f"Warning: {filepath} not found.")
            continue
        with open(filepath) as f:
            for line in f:
                target_id = line.strip()  #
                if target_id:
                    mapping_cono[target_id] = classification
    return mapping_cono


def cono_classification(tsv_file, species_id, mapping_dir="input/cono_hits"):
    # Look through Conotoxin Foldseek file and find matches for each node
    mapping_cono = cono_map(
        mapping_dir="input/classification_conotoxin",
        file_mapping=superfam_map,
    )

    # Prepare accumulator for heatmap
    species_list = species_id
    superfamilies = list(set(mapping_cono.values()))
    accumulator = {
        sf: {sp: {"total": 0.0, "count": 0} for sp in species_list}
        for sf in superfamilies
    }

    # Iterate through .tsv to build per-node classification
    cono_hits = {}
    df = pd.read_csv(
        tsv_file, sep="\t", header=None, names=["node", "target", "evalue"]
    )
    for _, row in df.iterrows():
        node, target, e_val = row["node"], row["target"], row["evalue"]
        if target in mapping_cono:
            sf = mapping_cono[target]
            entry = {"conotoxin_model": target, "classification": sf, "evalue": e_val}
            cono_hits.setdefault(node, []).append(entry)

            # Accumulate for heatmap
            sp = get_species(node, species_id)
            cell = accumulator[sf][sp]
            cell["total"] += e_val
            cell["count"] += 1

    # Build heatmap of average E-values
    cono_heatmap = {}
    for sf, sp_dict in accumulator.items():
        cono_heatmap[sf] = {}
        for sp, stats in sp_dict.items():
            if stats["count"] > 0:
                cono_heatmap[sf][sp] = stats["total"] / stats["count"]
            else:
                cono_heatmap[sf][sp] = None

    return cono_hits, cono_heatmap


def vzone_classification(file_path, species_id, mapping_dir="vzone_hits"):
    # Load the file mapping (3 columns: query, target, evalue)
    mapping = vzone_map(mapping_dir, toxinfold_map)

    # Retrive full species list
    species_list = species_id

    # Prepare accumulator with each cell holding {"total": 0, "count": 0}
    accumulator = {
        cl: {sp: {"total": 0.0, "count": 0} for sp in species_list}
        for cl in list(set(mapping.values()))
    }
    vzone_hits = {}
    df_af = pd.read_csv(
        file_path, sep="\t", header=None, names=["node", "target", "evalue"]
    )
    for _, row in df_af.iterrows():
        node = row["node"]
        target = row["target"]
        e_val = row["evalue"]
        # Look up the target in the mapping loaded from files
        if target in mapping:
            gen_class = mapping[target]
            entry = {"af_model": target, "classification": gen_class, "evalue": e_val}
            vzone_hits.setdefault(node, []).append(entry)

            # Get species from the node
            sp = get_species(node, species_id)

            if gen_class not in accumulator:
                accumulator[gen_class] = {
                    s: {"total": 0.0, "count": 0} for s in species_list
                }
            # Accumulate the E-value
            cell = accumulator[gen_class][sp]
            cell["total"] += e_val
            cell["count"] += 1
    # Compute the avg. E-value for each cell
    vzone_heatmap = {}
    for cl, species_dict in accumulator.items():
        vzone_heatmap[cl] = {}
        for sp, stats in species_dict.items():
            if stats["count"] > 0:
                vzone_heatmap[cl][sp] = (
                    stats["total"] / stats["count"]
                )  # average E-value
            else:
                vzone_heatmap[cl][sp] = None
    return vzone_hits, vzone_heatmap


# Visualization using cluster positions
# Make sure to list variables in the same order as in call to visualize_graph
def visualize_graph(
    G,
    query_targets,
    qt_pairs,
    output_json="data.json",
    vzone_hits=None,
    afdb_hits=None,
    cono_hits=None,
    cono_heatmap=None,
    vzone_heatmap=None,
    node_clusters=None,
):
    net = Network(notebook=True, width="1000px", height="800px", directed=False)

    # Predefine species colors for visualization - will prepend them as node colors
    species_colors = {
        "test_platymeris": "#882255",
        "test_corixa": "#AA4499",
        "test_vespa": "#AA3377",
        "test_myrmica": "#CC6677",
        "test_latrohesperus": "#661100",
        "test_euscorpius": "#C14444",
        "test_scolopendra": "#A52A2A",
        "test_xibalbanus": "#B03A48",
        "test_blenny": "#D46A6A",
        "test_sting": "#EE8866",
        "test_shrew": "#8C8C00",
        "test_crotalus": "#DDCC77",
        "test_najanaja": "#CCBB44",
        "test_heloderma": "#999933",
        "test_hapalochlaena": "#117733",
        "test_crassispira": "#228833",
        "test_unedogemmula": "#44AA99",
        "test_gemmula": "#66CCEE",
        "test_gloriamaris": "#88CCEE",
        "test_sponsalis": "#77AADD",
        "test_geographus": "#332288",
        "test_actinia": "#004488",
    }

    # Build dictionaries for PDB structures
    # Each node will have a PDB file in the AllSpecies directory
    query_pdb_dict = {}
    target_pdb_dict = {}
    for node in G.nodes():
        pdb_path = f"AllSpecies/{node}.pdb"
        if os.path.exists(pdb_path):
            with open(pdb_path, "r") as f:
                pdb_content = f.read()
            query_pdb_dict[node] = pdb_content
            target_pdb_dict[node] = pdb_content
        else:
            query_pdb_dict[node] = ""
            target_pdb_dict[node] = ""

    # Define a directory to store precomputed superpositions
    SUPERPOSED_DIR = "Superposed"
    if not os.path.exists(SUPERPOSED_DIR):
        os.makedirs(SUPERPOSED_DIR)

    # Precompute superposed PDBs for each query-target pair.
    query_superposed_pdb_paths = {}
    if not bare_mode:
        for key in qt_pairs:
            q, t = key.split("||")
            qpdb = query_pdb_dict.get(q, "")
            tpdb = query_pdb_dict.get(t, "")
            # Create a unique filename for this pair
            filename = f"{q}__{t}.pdb"  # Using __ as seperator
            file_path = os.path.join(SUPERPOSED_DIR, filename)

            if not os.path.exists(file_path):
                if qpdb and tpdb:
                    # Compute the superposition
                    aligned = superpose(qpdb, tpdb)
                    with open(file_path, "wt") as f_out:
                        f_out.write(aligned)
                else:
                    # If one of the PDBs is missing, write query pdb only
                    with open(file_path, "w") as f:
                        f.write(qpdb)

            # Store the relative path
            query_superposed_pdb_paths[key] = (
                f"/superposed/{os.path.basename(file_path)}"
            )
    else:  # If in bare mode, expect superposed files to be already present
        for key in qt_pairs:
            q, t = key.split("||")
            filename = f"{q}__{t}.pdb"
            file_path = os.path.join(SUPERPOSED_DIR, filename)
            if os.path.exists(file_path):
                query_superposed_pdb_paths[key] = (
                    f"/superposed/{os.path.basename(file_path)}"
                )
            else:
                query_superposed_pdb_paths[key] = query_pdb_dict.get(
                    q, ""
                )  # Fallback to query PDB if not found

    afdb_hits = load_AFDBhits("input/AFDBhits.tsv")

    # Prepare JSON data to be saved on disk
    # These data will be appended in the scripts.js file
    data = {
        "query_targets": query_targets,
        "qt_pairs": qt_pairs,
        "nodes": list(G.nodes()),
        "node_clusters": node_clusters,
        "query_superposed_pdb": query_superposed_pdb_paths,
        "vzone_hits": vzone_hits if vzone_hits else {},
        "vzone_heatmap": vzone_heatmap if vzone_heatmap else {},
        "afdb_hits": afdb_hits,
        "species_colors": species_colors,
        "cono_heatmap": cono_heatmap,
        "cono_hits": cono_hits,
    }

    with open(output_json, "w") as f:
        json.dump(data, f)

    # Group nodes by cluster
    clusters = {}
    for node in G.nodes():
        cl = G.nodes[node].get("cluster", 0)
        clusters.setdefault(cl, []).append(node)

    sorted_cluster_ids = sorted(clusters.keys())

    # Define grid parameters for clusters
    # Change setting to adjust the layout
    clusters_per_row = 10
    cluster_block_size = 500

    # Dictionary for block offsets
    cluster_offsets = {}
    for i, cl in enumerate(sorted_cluster_ids):
        col = i % clusters_per_row  # Compute column index
        row = i // clusters_per_row  # Compute row index
        offset_x = col * cluster_block_size
        offset_y = row * cluster_block_size
        cluster_offsets[cl] = (offset_x, offset_y)

    # Place nodes in each cluster using a grid layout
    node_spacing = 50

    for cl, nodes in clusters.items():
        num_nodes = len(nodes)
        # Compute number of columns required for grid layout
        grid_cols = math.ceil(math.sqrt(num_nodes))
        grid_rows = math.ceil(num_nodes / grid_cols)
        base_x, base_y = cluster_offsets.get(cl, (0, 0))

        for idx, node in enumerate(nodes):
            # Compute row and column for this node in the grid
            r = idx // grid_cols
            c = idx % grid_cols

            offset_node_x = c * node_spacing
            offset_node_y = r * node_spacing

            total_grid_width = (grid_cols - 1) * node_spacing
            total_grid_height = (grid_rows - 1) * node_spacing
            center_offset_x = (cluster_block_size - total_grid_width) / 2
            center_offset_y = (cluster_block_size - total_grid_height) / 2

            x = base_x + center_offset_x + offset_node_x
            y = base_y + center_offset_y + offset_node_y

            species = G.nodes[node]["species"]
            net.add_node(
                node,
                label=" ",
                color=species_colors.get(species, "#808080"),
                title=f"Protein: {node}",
                x=x,
                y=y,
                size=10,
                species=species,
            )

    # Add edges
    for edge in G.edges():
        net.add_edge(edge[0], edge[1])

    # Set options for the network visualization - disable physics for static layout (debugging only - no layout)
    net.set_options('{"physics": {"enabled": false}}')
    # Generate HTML for the network
    network_html = net.generate_html()

    with open("template.html") as f:
        html = f.read()
    species_colors_inject = (
        f"<script>const speciesColors = {json.dumps(species_colors)};</script>\n"
    )
    html = html.replace("<!-- SPECIES_COLORS_PLACEHOLDER -->", species_colors_inject)
    html = html.replace("<!-- CLUSTERING_PLACEHOLDER -->", network_html)
    with open("index.html", "w") as f:
        f.write(html)

    print(f"Interactive HTML graph saved as index.html and data as {output_json}")


def main(
    foldseek_files,
    species_id,
    bare_mode=False,
    chunk_size=100000,
):
    df = load_foldseek(foldseek_files, species_id=species_id, chunk_size=chunk_size)
    G = build_graph(df)
    components = list(nx.connected_components(G))
    node_clusters = {}
    for i, comp in enumerate(components):
        for node in comp:
            node_clusters[node] = i
    print(
        "Cluster counts (from connected components):",
        {i: len(comp) for i, comp in enumerate(components)},
    )
    print(f"Total clusters: {len(components)}")
    print(f"Total nodes: {sum(len(comp) for comp in components)}")

    # Build mapping from each query to its list of targets and pair metrics.
    query_targets = {}
    qt_pairs = {}
    for _, row in df.iterrows():
        query = row["query"]
        target = row["target"]
        query_targets.setdefault(query, []).append(target)
        key = f"{query}||{target}"
        qt_pairs[key] = {}

    # Load the AF classification foldseek file.
    vzone_hits_file = "input/vzonehits.tsv"
    vzone_hits = {}
    vzone_heatmap = {}
    if os.path.exists(vzone_hits_file):
        vzone_hits, vzone_heatmap = vzone_classification(
            vzone_hits_file,
            species_id,
            mapping_dir="input/classification_vzone/",
        )
    cono_hits_file = "input/conotoxinhits.tsv"
    cono_hits = {}
    cono_heatmap = {}
    if os.path.exists(cono_hits_file):
        cono_hits, cono_heatmap = cono_classification(
            cono_hits_file,
            species_id,
            mapping_dir="input/classification_conotoxin/",
        )

    visualize_graph(
        G,
        query_targets,
        qt_pairs,
        "data.json",
        vzone_hits=vzone_hits,
        cono_hits=cono_hits,
        cono_heatmap=cono_heatmap,
        vzone_heatmap=vzone_heatmap,
        node_clusters=node_clusters,
    )
    with open(os.path.join(os.path.dirname(__file__), "node_clusters.json"), "w") as f:
        json.dump(node_clusters, f)
    print("Interactive HTML graph saved as index.html with data in data.json")


# This script is intended to be run from the command line
# Usage: python CoNAT.py --bare, if you want to skip superposition step
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Network visualization tool")
    parser.add_argument(
        "--bare",
        action="store_true",
        help="Run in barebones mode (skip superposition step - faster, but no superposed PDBs)",
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=5000,
        help="Number of rows per chunk to process",
    )
    args = parser.parse_args()
    bare_mode = args.bare

    # Needs to be updated according to ColabFold file identifiers - renaming was used for these.
    species_id = [
        "test_crotalus",
        "test_sponsalis",
        "test_sting",
        "test_vespa",
        "test_shrew",
        "test_latrohesperus",
        "test_najanaja",
        "test_heloderma",
        "test_myrmica",
        "test_geographus",
        "test_blenny",
        "test_euscorpius",
        "test_gloriamaris",
        "test_platymeris",
        "test_scolopendra",
        "test_corixa",
        "test_unedogemmula",
        "test_gemmula",
        "test_crassispira",
        "test_xibalbanus",
        "test_actinia",
        "test_hapalochlaena",
    ]
    foldseek_files = glob.glob("input/result.tsv")
    classification_file = "input/afdbhits.tsv"

    main(
        foldseek_files,
        species_id,
        bare_mode=bare_mode,
        chunk_size=args.chunk_size,
    )
