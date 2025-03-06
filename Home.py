import streamlit as st
import pandas as pd
import plotly.express as px

# Page Configuration
st.set_page_config(page_title="Explore", layout='wide')
st.title("NuoHMMER Web v1b")

# Overview Section
with st.expander("Info", expanded=True):
    st.markdown("""
    **NuoHMMER Web** is an interactive platform to explore the distribution and diversity of Complex I subunits 
    across prokaryotic genomes. The dataset includes over 45,000 genomes, analyzing Complex I subunits presence and its variations.
    """)

# Load Data (Optimized with Caching)
@st.cache_data
def load_data():
    data_file = "https://raw.githubusercontent.com/akshayonly/demostreamlit/refs/heads/main/data/data.csv"
    return pd.read_csv(data_file)

# Tabs
tab1, tab2 = st.tabs(["Overview", "Explore"])

# Colors for the Pie Chart
colors = ["#78C679", "#FEC44F", "#FEE08B", "#57A4E3", "#D0D1E6", "#E41A1C"]

# **Tab 1: Overview**
with tab1:
    df = pd.DataFrame({
        "Nuo subunit status": [
            "Complete C-I", "Complete C-I (CD fused)", "Complete C-I (BCD fused)",
            "C-I like (EF/EFG missing)", "Incomplete C-I", "Nuo subunits absent"
        ],
        "Species": [3884, 1363, 11, 415, 2171, 3211]
    })

    fig = px.pie(df, names="Nuo subunit status", values="Species", title="Species Distribution",
                 color="Nuo subunit status", color_discrete_sequence=colors)
    fig.update_traces(textposition='outside', textinfo='label+percent')
    fig.update_layout(showlegend=False)

    st.plotly_chart(fig)

# **Tab 2: Explore**
with tab2:
    data = load_data()

    # Column Selection
    subunits = ['NuoA', 'NuoB', 'NuoBCD', 'NuoC', 'NuoCD', 'NuoD', 'NuoE', 'NuoF', 
                'NuoG', 'NuoH', 'NuoI', 'NuoJ', 'NuoK', 'NuoL', 'NuoM', 'NuoN']
    
    selected_columns = ['Assembly Accession', 'Reference', 'Accession', 'Replicon', 'Variation', 'Count'] + subunits

    # Species Selection
    selected_value = st.selectbox("Choose Species:", sorted(data['Species'].dropna().astype(str).unique()))

    # Filtering Data
    filtered_data = data[data['Species'] == selected_value].reset_index(drop=True)

    if filtered_data.empty:
        st.warning(f"No data available for '{selected_value}'.")
    else:
        # Extract Taxonomy Information
        taxonomy_info = filtered_data[["Superkingdom", "Phylum", "Class", "Family"]].iloc[0].to_dict()

        col1, col2, col3 = st.columns(3)
        for col, (label, value) in zip([col1, col2, col3], taxonomy_info.items()):
            col.markdown(f"""
                <div style="font-size:14px; font-weight:normal; text-align:left;">
                    {label}<br><span style="font-size:16px; font-weight:bold; color:#27ae60;">{value}</span>
                </div>
            """, unsafe_allow_html=True)

        st.markdown("---")

        # Complex I Analysis
        total_genomes = filtered_data['Assembly Accession'].nunique()
        species_with_complex_i = filtered_data[filtered_data['Variation'] != 'No Complex I']['Assembly Accession'].nunique()
        species_without_complex_i = total_genomes - species_with_complex_i

        col1, col2, col3 = st.columns(3)
        col1.metric(label="Genomes Records", value=f"{total_genomes}")
        col2.metric(label="Genomes with Complex I subunit(s)", value=f"{species_with_complex_i}")
        col3.metric(label="Genomes without Complex I subunit(s)", value=f"{species_without_complex_i}")

        st.markdown("---")

        # Check for Functioning Complex I
        species_variations = set(filtered_data['Variation'].unique())
        functioning_variations = {'Nuo12', 'Nuo13', 'Nuo14', 'Nuo14-EF', 'Nuo14-EFG'}
        present_variations = functioning_variations & species_variations

        functioning_ci = "Present" if present_variations else "Absent"

        col1, col2, col3 = st.columns(3)
        col1.metric(label="Functioning Complex I", value=functioning_ci)
        col2.metric(label="Complex I Variation", value="; ".join(present_variations) if present_variations else "None")
        col3.metric(label="Total Subunits", value=filtered_data[subunits].sum().ge(1).sum())

        st.markdown("---")

        # Display Data Table
        if st.checkbox(f'Display {total_genomes} records'):
            st.dataframe(filtered_data[selected_columns], use_container_width=True)

# Footer
st.markdown("---")
st.caption("Developed at TIFR Mumbai by [Anand Research Group](https://www.anandresearch.com/) ", unsafe_allow_html=True)
