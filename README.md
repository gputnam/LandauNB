# Jupyter Notebook for Running Diffuse MPV 

Installation instructions:

Pre-Requisites: ROOT installation with MathMore enabled. Version of
python3.

python -m venv env
source env/bin/activate
pip install -r requirements.txt
pip install -r requirelandau.txt
python -m ipykernel install --name env --user

# Run the Jupyter notebook
jupyter notebook # open the DiffuseMPV.ipynb notebook and set the Kernel to env
