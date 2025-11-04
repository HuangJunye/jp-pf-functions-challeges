# QDC 2025 Challenge Installation Guide

## Step 0: Create your IBM Cloud account.
If you haven't done so already, set up your IBM Cloud account by following the instructions [here](https://quantum.cloud.ibm.com/docs/en/guides/cloud-setup).

## Step 1: Clone repository.

In your terminal, navigate to a folder where you would like to place the QDC 2025 Challenge repository:
```
cd </path/to/your/folder>
```
Then, clone the repository by running:
```
git clone https://github.com/qiskit-community/qdc-challenges-2025.git
```

## Step 2: Create virtual environment.

Please use Python 3.11-3.13. If you do not have one of these versions of Python installed on your machine, please follow the instructions [here](https://www.python.org/downloads/). Then, run the following code in your open terminal:
```
python3 -m venv qdc2025-venv
source qdc2025-venv/bin/activate
```
to create and activate a new virtual environment named **`qdc2025-venv`**.

## Step 3: Install required packages.

Finally, install the required packages by running the following line in your open terminal:
```
cd qdc-challenges-2025
pip install -r requirements.txt
```
Note that some challenges use `graphviz` for plotting, which needs to be installed on your machine independently of Python environment. If you do not have `graphviz` on your machine, please follow instructions [here](https://graphviz.org/download/). 

## Step 4: Register the virtual environment as a Jupyter kernel.

If you plan to work in JupyterLab, run the following lines to register your virtual environment as a Jupyter kernel and launch Jupyter. Once Jupyter launches and you open a challenge notebook, please click **`Kernel > Change Kernel...`** and then select **`QDC 2025`** from the dropdown menu to use the correct virtual environment.
```
python -m ipykernel install --user --name=qdc2025-venv --display-name "QDC 2025"
jupyter notebook
```

## Step 5: Save your token. 

Please navigate to `save_account.ipynb` and run the notebook to save your token. This must be done correctly in order for your challenge submissions to be graded.

You're ready to get started! Good luck :)
