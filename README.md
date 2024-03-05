# Guide to Using PMC for Polymer Property Prediction

Welcome! This guide will help you use the Persistent Multi-cover (PMC) machine learning approach for predicting polymer properties. We'll walk you through obtaining the data, setting up your environment, and executing the code to generate PMC features from your dataset.

## Data Availability

You can find all the datasets used in this project in the `/data/` directory. This includes both the original and augmented datasets.

## Getting Started

Before diving into the PMC feature generation, ensure you have the following Python packages installed:

- `numpy`
- `pandas`
- `sentence_transformers`
- `os`
- `re`
- `copy`
- `math`
- `matplotlib`
- `seaborn`
- `rdkit` (Version = 2023.03.3)
- `gudhi` (Version = 3.8.0)

Additionally, you'll need the `rhomboidtiling` software, available at https://github.com/geoo89/rhomboidtiling

## How to Use

Follow these steps to prepare and process your dataset:

1. **Prepare Your Dataset**: Place your dataset's CSV file in the `/data` directory. 

2. **Configure Settings**: Update the `config.yaml` file with your specific details. This includes:
   - `dataset_name`: The name of your dataset file (omit the '.csv' extension).
   - `save_path_pre`: The location where you want to save the coordinates and filtration of Delaunay slice file.
   - `rhomboidtiling_path`: The installation directory of the rhomboidtiling software. Ensure this directory contains the 'main' file, indicating a correct installation.

3. **Running the Code**:
   - Without data augmentation: Run **Multicover_process1.py** followed by **Multicover_process2.py**.
   - With data augmentation: Run **Multicover_aug_process1.py** followed by **Multicover_aug_process2.py**.

By following these steps, you'll be set to generate PMC features for your dataset in the '/result/' directory.

## Visualization

We also provide code snippets for visualizing different aspects of the PMC process, such as Multi-cover, Delaunay slice, and Rhomboid tiling. Here’s how to use them:

### Multi-cover Visualization

1. **Multi-cover (2D Visualization)**: To visualize Multi-cover in two dimensions, run **Multi-cover_visualization.py**. Remember to adjust the 'Vertexset' value within the code to match the coordinates of your point set.

### Rhomboid Tiling Visualization

2. **Rhomboid Tiling**:
   - To visualize Rhomboid tiling, you'll need to use the `rhomboidtiling` software with the 'cbifi' and 'rhomboids' options to generate two more file from your oringinal point set file first, then run **rhomboidtiling_visualization.py**
   - As an example:
     - `4points.txt` is the original coordinates of the point set file.
     - `4points_rhomboids.txt` is generated using the 'rhomboids' option.
     - `4points_rhomboidtiling.txt` is generated using the 'cbifi' option.
