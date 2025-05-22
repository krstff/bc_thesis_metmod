# MetMod: A Workflow for GEM Reconstruction

MetMod is a lightweight Python library and workflow designed to support **genome-scale metabolic model (GEM)** reconstruction, especially for **non-model organisms**. This repository contains the core implementation, documentation, and interactive notebooks demonstrating its practical use.

## About

Genome-scale metabolic models provide a detailed map of an organism's metabolism and are widely used in systems biology and biotechnology. However, reconstructing these models—particularly for non-model organisms—can be difficult due to limited data, annotation inconsistencies, and integration challenges.

This repository accompanies the [Bachelor's thesis]() _"Computational Workflow for Genome-Scale Metabolic Model Reconstruction of Non-Model Organisms"_ and presents a functional reconstruction workflow addressing those challenges. It includes:

- Tools for data integration and consistency
- Visualization support (via Escher)
- Support for model evaluation and refinement
- Practical examples on a real draft model

## Repository Contents

- `metmod.py`: The source code of the MetMod library.
- `workflow.ipynb`: A Jupyter notebook showcasing key MetMod functions.
- `workflow_draft.ipynb`: Workflow demonstration on a draft model.
- `workflow_core.ipynb`: Workflow demonstration on a curated core model.
- `docs_metmod/`: Documentation directory for MetMod.
- `inputs/`: Input models.
- `models/`: Finalized/polished metabolic models.
- `maps/`: Escher map files used for pathway visualization, downloaded from https://escher.github.io/.
- `expression_data/`: Gene expression data used in the workflows.
- `additional/`: Additional files, including a Cytoscape session and exported maps.

## Non-Public Data

Due to publication restrictions, some data used in the case study is hosted in a [private repository](https://gitlab.fi.muni.cz/xmatust/bc_thesis_non_public_data):

- `caldimonas_core_model.sbml` - A manually curated core model by Rosputinský
- `caldimonas_core_map.json` - An Escher map created from the model
- `RNASeq2_RAW_CountTable.csv` - RNA-seq data measured by Musilová

## Getting Started

To run the workflow:

```bash
# Create a new environment
conda create --name metmod-test python=3.9 -y
conda activate metmod-test
conda install pip -y
# Install Jupyter into the environment and install dependencies
conda install jupyter
pip install escher notebook
jupyter labextension install @jupyter-widgets/jupyterlab-manager
jupyter labextension install escher
# Install Python requirements
pip install -r requirements.txt
# Run the .ipynb notebooks
