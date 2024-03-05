# Dataset Overview

This dataset comprises a comprehensive collection of monomer data, meticulously organized into three columns. Each row represents a unique monomer, providing essential details for researchers and practitioners working in the field of chemistry, materials science, and related disciplines. Below is an overview of the dataset structure:

## Columns Description:

**smiles**: The SMILES string of the monomer.

**index**: A unique identifier assigned to each polymer associated with the monomer within the dataset.

**value**: The label associated with each monomer. 

We also provide the dataset after augmentation. The files contain four column:

## Columns Description:

**smiles**: The SMILES string of the monomer after augmentation.

**index**: A unique identifier assigned to each polymer associated with the monomer within the dataset.

**value**: The label associated with each monomer. 

**processed_strings**: The SMILES string of the monomer where the placeholder * has been meticulously replaced with the actual atom that occupies the subsequent position in the polymer.
