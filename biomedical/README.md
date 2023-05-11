# Brain Tumour Classification Project

This project aims to classify brain tumour types based on histopathological data. We utilize the "Digital Brain Tumour Atlas" dataset, which contains information about a variety of tumour types. For now, the project includes preprocessing the data, performing exploratory data analysis, and implementing a convolutional neural network for tumour type classification.

## Installation

This project used Conda to manage the environment. Get it [here](https://www.anaconda.com/download).

You can create a new conda environment with the following command:

```bash
conda env create -f environment.yml
```

Activate using 

```bash
conda activate brain-tumour-classification
```

## Project Structure
Here's a brief overview of the important files/folders:

`annotation.csv`: This CSV file contains metadata about the patients and tumors.

`Data preprocessing.ipynb`: This Jupyter notebook is used for loading the annotation data, defining a script for downloading the histopathological data, and preprocessing it.

`Data analysis.ipynb`: This Jupyter notebook is used for performing an exploratory data analysis on the preprocessed data, and training a convolutional neural network to classify tumor types.

`processed/`: This directory contains the processed histopathological data.

## Download data

To download and process the data, you need to get an authorization token. To get it, first [request access to the dataset](https://data-proxy.ebrains.eu/datasets/8fc108ab-e2b4-406f-8999-60269dc1f994). Then inspect the network traffic in the browser (F12 for Chrome) and search for the header `Authorization`. You can then insert the token value into the corresponding variable in `Data preprocessing.ipynb`.