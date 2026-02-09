import streamlit as st
import pandas as pd
import altair as alt
import requests
from stmol import showmol
import py3Dmol

# --- 1. CONFIGURATION ---
st.set_page_config(layout="wide", page_title="Bio-Sequencer Pro", page_icon="🧬")

# Custom CSS for "Billionaire Professional" Look
st.markdown("""
<style>
    .stTextArea textarea {
        font-family: 'Courier New', monospace;
        background-color: #f0f2f6;
        color: #31333F;
    }
    .step-title {
        font-size: 1.2rem;
        font-weight: bold;
        color: #00C896; /* Teal Accent */
        margin-bottom: 10px;
    }
</style>
""", unsafe_allow_html=True)

# --- 2. HEADER ---
st.title("🧬 Genomic Sequence Analyzer")
st.markdown("""
**From Code to Structure:** Transform raw DNA sequences into biological insights and 3D protein structures using Computational Biology.
""")

# --- 3. SIDEBAR (INPUT) ---
st.sidebar.header("1. Input Data")

default_seq = ">Example_Sequence\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
input_mode = st.sidebar.radio("Input Mode", ["Paste Sequence", "Use Example"])

if input_mode == "Use Example":
    raw_sequence = default_seq
else:
    raw_sequence = st.sidebar.text_area("Paste DNA Sequence (FASTA format supported)", height=250, placeholder=">Sequence_1\nATGC...")

# --- 4. LOGIC ENGINE ---

def clean_sequence(seq):
    """Removes FASTA headers and whitespace"""
    lines = seq.splitlines()
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
    return "".join(lines).upper().replace(" ", "")

def transcribe_dna(dna):
    """Replaces T with U"""
    return dna.replace("T", "U")

def translate_rna(rna):
    """Translates mRNA to Amino Acids using Standard Genetic Code"""
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = ""
    for i in range(0, len(rna), 3):
        if i+3 <= len(rna):
            codon = rna[i:i+3]
            protein += codon_table.get(codon, 'X')
    return protein

def get_3d_structure(protein_sequence):
    """Fetches PDB data from ESMFold API (Meta AI)"""
    # API Endpoint for Protein Folding
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    response = requests.post(url, data=protein_sequence, verify=False)
    if response.status_code == 200:
        return response.text
    else:
        return None

# PROCESS DATA
cleaned_dna = clean_sequence(raw_sequence)
mrna_seq = transcribe_dna(cleaned_dna)
protein_seq = translate_rna(mrna_seq)

# --- 5. MAIN DISPLAY (Tabs) ---

tab1, tab2, tab3, tab4 = st.tabs(["1. Quality Check", "2. Transcription", "3. Translation", "4. 3D Structure"])

with tab1:
    st.markdown('<div class="step-title">Step 1: Sequence Cleaning</div>', unsafe_allow_html=True)
    col1, col2 = st.columns(2)
    with col1:
        st.caption("Raw Input")
        st.code(raw_sequence, language="text")
    with col2:
        st.caption("Cleaned DNA")
        st.code(cleaned_dna, language="text")
        
    # Nucleotide Counts Chart
    counts = {'A': cleaned_dna.count('A'), 'T': cleaned_dna.count('T'), 'G': cleaned_dna.count('G'), 'C': cleaned_dna.count('C')}
    df = pd.DataFrame.from_dict(counts, orient='index', columns=['Count']).reset_index()
    chart = alt.Chart(df).mark_bar().encode(
        x='index', y='Count', color=alt.value("#00C896")
    ).properties(height=300)
    st.altair_chart(chart, use_container_width=True)

with tab2:
    st.markdown('<div class="step-title">Step 2: Transcription (DNA → mRNA)</div>', unsafe_allow_html=True)
    st.write("Thymine (T) is replaced by Uracil (U).")
    st.code(mrna_seq, language="text")

with tab3:
    st.markdown('<div class="step-title">Step 3: Translation (mRNA → Protein)</div>', unsafe_allow_html=True)
    st.write(f"Total Amino Acids: {len(protein_seq)}")
    st.code(protein_seq, language="text")

with tab4:
    st.markdown('<div class="step-title">Step 4: AI Folding Prediction</div>', unsafe_allow_html=True)
    st.write("Using **ESMFold (Meta AI)** to predict 3D structure.")
    
    if st.button("🧬 Generate 3D Structure"):
        with st.spinner("Contacting Meta AI... Folding Protein..."):
            pdb_data = get_3d_structure(protein_seq)
            
        if pdb_data:
            st.success("Structure Generated!")
            view = py3Dmol.view(width=800, height=500)
            view.addModel(pdb_data, "pdb")
            view.setStyle({'cartoon': {'color': 'spectrum'}})
            view.zoomTo()
            showmol(view, height=500, width=800)
        else:
            st.error("Error fetching structure. API might be busy.")
    else:
        st.info("Click the button above to start.")