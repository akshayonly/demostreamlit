import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path

# Set the page layout to wide
st.title("NuoHMMER Web v1b")
st.info('This site is actively being improved, and new features may be added. Core functionalities are stable for research use.')
st.set_page_config(layout="wide")

# Corrected file paths using GitHub Raw URLs
data_file = "https://raw.githubusercontent.com/akshayonly/demostreamlit/main/data/data.csv"
metadata_file = "https://raw.githubusercontent.com/akshayonly/demostreamlit/main/data/metadata.csv"

# Check if required files exist (Skip local path existence check since we're using URLs)
try:
    data = pd.read_csv(data_file, encoding='utf-8')
    metadata = pd.read_csv(metadata_file, encoding='utf-8')
except Exception as e:
    st.error(f"Error loading data files: {e}")
    st.stop()

# Search Options using Radio Button
search_option = st.radio(
    'Choose a search criteria:',
    ['NCBI Accession Number', 'Species']
)

# Input and filtering based on search criteria
selected_value = None
filtered_data = pd.DataFrame()

if search_option == 'NCBI Accession Number':
    selected_value = st.selectbox('Choose an Accession Number or type to search', sorted(data['Accession'].unique()))
    filtered_data = data[data['Accession'] == selected_value]

else:  # search_option == 'Species'
    selected_value = st.selectbox('Choose a Species or type to search', sorted(data['Species'].unique()))
    filtered_data = data[data['Species'] == selected_value]

# Display a message if no data is available for the selected value
if filtered_data.empty:
    st.text(f"No data available for the selected '{selected_value}' under '{search_option}'.")
else:
    number_species = filtered_data['Species'].nunique()
    species_with_complex_i = filtered_data[filtered_data['Variation'] != 'Absent']
    species_without_complex_i = filtered_data[(filtered_data['Variation'] == 'Absent') &
                                              (~filtered_data['Species'].isin(species_with_complex_i['Species']))]

    num_with_complex_i = species_with_complex_i['Species'].nunique()
    num_without_complex_i = species_without_complex_i['Species'].nunique()

    # Display human-friendly results
    if num_with_complex_i == number_species:
        st.text(f"All {number_species} species under '{selected_value}' have Complex I Subunit(s).")
    elif num_without_complex_i == number_species:
        st.text(f"All {number_species} species under '{selected_value}' lack Complex I Subunit(s).")
    else:
        st.text(f"Under '{selected_value}', there are {number_species} species in total.")
        st.text(f"{num_with_complex_i} species have Complex I Subunit(s), while {num_without_complex_i} species do not.")

    # Display filtered data
    st.subheader('Species with Complex I Subunit(s)')
    st.dataframe(species_with_complex_i, use_container_width=True)

    st.subheader('Species without Complex I Subunit(s)')
    st.dataframe(species_without_complex_i, use_container_width=True)
