import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle
import glob
from padelpy import padeldescriptor

def pubchem():
    xml_files = glob.glob("*.xml")  
    xml_files.sort()
    xml_files
    FP_list = ['AtomPairs2DCount',
    'AtomPairs2D',
    'EState',
    'CDKextended',
    'CDK',
    'CDKgraphonly',
    'KlekotaRothCount',
    'KlekotaRoth',
    'MACCS',
    'PubChem',
    'SubstructureCount',
    'Substructure']
    fp = dict(zip(FP_list, xml_files))
    fp['PubChem']
    fingerprint = 'PubChem'
    fingerprint_output_file = ''.join([fingerprint,'.csv']) #Substructure.csv
    fingerprint_descriptortypes = fp[fingerprint]

# Molecular descriptor calculator
def desc_calc():    
    xml_files = glob.glob("*.xml")  
    xml_files.sort()
    FP_list = ['AtomPairs2DCount',
    'AtomPairs2D',
    'EState',
    'CDKextended',
    'CDK',
    'CDKgraphonly',
    'KlekotaRothCount',
    'KlekotaRoth',
    'MACCS',
    'PubChem',
    'SubstructureCount',
    'Substructure']
    fp = dict(zip(FP_list, xml_files))

    fp['PubChem']
    fingerprint = 'PubChem'
    fingerprint = 'PubChem'
    fingerprint_output_file = ''.join([fingerprint,'.csv']) #Substructure.csv
    fingerprint_descriptortypes = fp[fingerprint]

    padeldescriptor(mol_dir='molecule.smi', 
                    d_file=fingerprint_output_file, #'Substructure.csv'
                    #descriptortypes='SubstructureFingerprint.xml', 
                    descriptortypes= fingerprint_descriptortypes,
                    detectaromaticity=True,
                    standardizenitro=True,
                    standardizetautomers=True,
                    threads=2,
                    removesalt=True,
                    log=True,
                    fingerprints=True) 


# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building
def build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('sert_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)

# Logo image
image = Image.open('logo.png')

st.image(image, use_column_width=True)

# Page title
st.markdown("""
# Bioactivity Prediction App (Serotonin)
This app allows you to predict the bioactivity towards inhibting the `Serotonin` enzyme. 
**Credits**
- App built in `Python` + `Streamlit` by [Chanin Nantasenamat](https://medium.com/@chanin.nantasenamat) (aka [Data Professor](http://youtube.com/dataprofessor))
- Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Read the Paper]](https://doi.org/10.1002/jcc.21707).
---
""")

# Sidebar
with st.sidebar.header('1. Upload your CSV data'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    st.sidebar.markdown("""
[Example input file](https://raw.githubusercontent.com/dataprofessor/bioactivity-prediction-app/main/example_acetylcholinesterase.txt)
""")

if st.sidebar.button('Predict'):
    load_data = pd.read_table(uploaded_file, sep=' ', header=None)
    load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)

    st.header('**Original input data**')
    st.write(load_data)

    with st.spinner("Calculating descriptors..."):
        desc_calc()

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors**')
    desc = pd.read_csv('PubChem.csv')
    st.write(desc)
    st.write(desc.shape)

    # Read descriptor list used in previously built model
    st.header('**Subset of descriptors from previously built models**')
    Xlist = list(pd.read_csv('descriptor_list.csv').columns)
    desc_subset = desc[Xlist]
    st.write(desc_subset)
    st.write(desc_subset.shape)

    # Apply trained model to make prediction on query compounds
    build_model(desc_subset)
else:
    st.info('Upload input data in the sidebar to start!')