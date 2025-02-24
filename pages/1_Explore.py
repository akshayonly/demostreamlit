import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path

# Set the page layout to wide
st.set_page_config(layout="wide")
st.title("NuoHMMER Web v1b")
st.info('This site is actively being improved, and new features may be added. Core functionalities are stable for research use.')

# Corrected file paths using GitHub Raw URLs
data_file = "https://raw.githubusercontent.com/akshayonly/demostreamlit/main/data/data.csv"
metadata_file = "https://raw.githubusercontent.com/akshayonly/demostreamlit/main/data/metadata.csv"
# seqs_file = "https://raw.githubusercontent.com/akshayonly/demostreamlit/main/data/nuo_searched_seqs.faa"

# Check if required files exist (Skip local path existence check since we're using URLs)
try:
    data = pd.read_csv(data_file, encoding='utf-8')
    metadata = pd.read_csv(metadata_file, encoding='utf-8')
except Exception as e:
    st.error(f"Error loading data files: {e}")
    st.stop()

# Title of the web app
st.title('NuoHMMER Results')

# Search Options using Radio Button
search_option = st.radio(
    'Choose a search criteria:',
    ['Genus', 'NCBI Accession Number', 'Species']
)

# Input and filtering based on search criteria
selected_value = None
filtered_data = pd.DataFrame()

if search_option == 'Genus':
    selected_value = st.selectbox('Choose a Genus or type to search', sorted(data['Genus'].unique()))
    filtered_data = data[data['Genus'] == selected_value]

elif search_option == 'NCBI Accession Number':
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

    # Display Subunit Sequences if checkbox is selected
    # if search_option == 'NCBI Accession Number' and st.checkbox('Display Subunit(s) Sequences'):
    #     sequence_data = metadata[metadata['Accession'] == selected_value][['ProteinAccession', 'Subunit']]
    #     prot_accessions = sequence_data.set_index('ProteinAccession')['Subunit'].to_dict()

    #     # Read sequence file from URL
    #     try:
    #         import requests
    #         seqs_response = requests.get(seqs_file)
    #         seqs_content = seqs_response.text

    #         combined_output = []
    #         for record in SeqIO.parse(seqs_content.splitlines(), 'fasta'):
    #             if record.id in prot_accessions:
    #                 cleaned_sequence = str(record.seq).replace('*', '')
    #                 seq_record = SeqRecord(
    #                     seq=cleaned_sequence,
    #                     id=record.id,
    #                     description=record.description
    #                 )
    #                 combined_output.append(f">{prot_accessions[record.id]} {seq_record.description}\n{seq_record.seq}")

    #         if combined_output:
    #             st.subheader('Subunit Sequences')
    #             st.text("\n".join(combined_output))
    #     except Exception as e:
    #         st.error(f"Error loading sequence file: {e}")
