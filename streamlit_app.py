import altair as alt
import pandas as pd
import py3Dmol
import requests
import streamlit as st
from stmol import showmol

st.set_page_config(layout="wide", page_title="Bio-Sequencer Pro", page_icon="🧬")

st.markdown("""
<style>
    :root { --ink: #17212b; --muted: #52616b; --teal: #008f83; --coral: #e76f51; --mist: #eef5f3; }
    .stApp, [data-testid="stAppViewContainer"] { background: #f7f9f8; color: var(--ink); }
    [data-testid="stHeader"] { background: transparent; }
    [data-testid="stMainBlockContainer"] { padding-top: 3.5rem; }
    [data-testid="stSidebar"] { background: #17212b; }
    [data-testid="stSidebar"] * { color: #f7f9f8; }
    [data-testid="stSidebar"] .stTextArea textarea { background: #253542; color: #f7f9f8; }
    .stTextArea textarea, code { font-family: 'Courier New', monospace; }
    .hero { border-left: 6px solid var(--coral); padding: 0.25rem 1.25rem; margin-bottom: 1.5rem; }
    .hero p { color: var(--muted); font-size: 1.05rem; }
    .step-title { color: var(--teal); font-size: 1.25rem; font-weight: 700; margin-bottom: 0.75rem; }
    [data-baseweb="tab-list"] { gap: 0.5rem; border-bottom: 1px solid #cbd8d4; }
    button[role="tab"] { color: var(--muted) !important; font-weight: 600; }
    button[role="tab"][aria-selected="true"] { color: var(--coral) !important; }
    [data-testid="stAlert"] { color: var(--ink); }
    [data-testid="stMetric"] { background: white; border: 1px solid #d7e2df; border-radius: 8px; padding: 0.75rem; }
    [data-testid="stMetricLabel"] { color: var(--muted); }
    .metric-strip { background: var(--mist); border-radius: 8px; padding: 0.75rem 1rem; }
    div.stButton > button { border-radius: 6px; border: 0; background: var(--teal); color: white; font-weight: 700; }
</style>
""", unsafe_allow_html=True)

CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

DEFAULT_SEQUENCE = ">Example_Sequence\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"


def clean_sequence(sequence):
    """Remove FASTA headers, whitespace, and line breaks."""
    return "".join(
        line.strip().upper()
        for line in sequence.splitlines()
        if not line.strip().startswith(">")
    )


def validate_sequence(sequence):
    errors = []
    if not sequence:
        errors.append("Enter a DNA sequence or choose the example sequence.")
    invalid = sorted(set(sequence) - set("ACGT"))
    if invalid:
        errors.append(f"Invalid DNA character(s): {', '.join(invalid)}.")
    if sequence and len(sequence) < 3:
        errors.append("The sequence must contain at least one complete codon.")
    if sequence and len(sequence) % 3:
        errors.append("The sequence length must be divisible by 3 for this reading frame.")
    return errors


def transcribe_dna(dna):
    return dna.replace("T", "U")


def translate_rna(rna):
    protein = []
    for index in range(0, len(rna) - 2, 3):
        amino_acid = CODON_TABLE[rna[index:index + 3]]
        if amino_acid == "*":
            break
        protein.append(amino_acid)
    return "".join(protein)


def get_3d_structure(protein_sequence):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        response = requests.post(
            url,
            data=protein_sequence,
            headers={"Content-Type": "text/plain"},
            timeout=(10, 120),
        )
        response.raise_for_status()
        return response.text, None
    except requests.exceptions.Timeout:
        return None, "The folding service took too long to respond. Try again with a shorter protein."
    except requests.exceptions.RequestException as error:
        return None, f"The folding service returned an error: {error}"


for key, value in {
    "active_raw": "",
    "cleaned_dna": "",
    "quality_passed": False,
    "transcription_passed": False,
    "translation_passed": False,
    "pdb_data": None,
    "fold_error": None,
}.items():
    st.session_state.setdefault(key, value)

st.markdown("""
<div class="hero">
    <h1>🧬 Genomic Sequence Analyzer</h1>
    <p>Move from a validated DNA sequence to mRNA, protein, and a predicted 3D structure.</p>
</div>
""", unsafe_allow_html=True)

st.sidebar.header("Start an analysis")
with st.sidebar.form("sequence_form"):
    input_mode = st.radio("Input source", ["Use Example", "Paste Sequence"], index=0)
    if input_mode == "Use Example":
        raw_input = DEFAULT_SEQUENCE
        st.caption("A validated teaching sequence is ready to run.")
    else:
        raw_input = st.text_area(
            "DNA sequence (FASTA supported)",
            height=220,
            placeholder=">Sequence_1\nATGGCC...",
        )
    analyze = st.form_submit_button("Analyze sequence", use_container_width=True)

if analyze:
    cleaned = clean_sequence(raw_input)
    errors = validate_sequence(cleaned)
    st.session_state.update({
        "active_raw": raw_input,
        "cleaned_dna": cleaned,
        "quality_passed": not errors,
        "transcription_passed": False,
        "translation_passed": False,
        "pdb_data": None,
        "fold_error": None,
    })
    if errors:
        for error in errors:
            st.sidebar.error(error)
    else:
        st.sidebar.success("Sequence validated. Continue through the stages.")

cleaned_dna = st.session_state.cleaned_dna
mrna_seq = transcribe_dna(cleaned_dna) if st.session_state.quality_passed else ""
protein_seq = translate_rna(mrna_seq) if st.session_state.transcription_passed else ""

st.sidebar.divider()
st.sidebar.caption("Workflow status")
status = [
    ("Quality Check", st.session_state.quality_passed),
    ("Transcription", st.session_state.transcription_passed),
    ("Translation", st.session_state.translation_passed),
    ("3D Structure", bool(st.session_state.pdb_data)),
]
for label, complete in status:
    st.sidebar.write(f"{'✓' if complete else '○'} {label}")

tab1, tab2, tab3, tab4 = st.tabs([
    "1. Quality Check", "2. Transcription", "3. Translation", "4. 3D Structure"
])

with tab1:
    st.markdown('<div class="step-title">Step 1 · Validate the DNA input</div>', unsafe_allow_html=True)
    if not st.session_state.active_raw:
        st.info("Choose an input source in the sidebar, then select Analyze sequence.")
    else:
        col1, col2 = st.columns(2)
        with col1:
            st.caption("Raw input")
            st.code(st.session_state.active_raw, language="text")
        with col2:
            st.caption("Cleaned DNA")
            st.code(cleaned_dna, language="text")
        counts = {base: cleaned_dna.count(base) for base in "ATGC"}
        gc_content = ((counts["G"] + counts["C"]) / len(cleaned_dna) * 100) if cleaned_dna else 0
        metric_col1, metric_col2, metric_col3 = st.columns(3)
        metric_col1.metric("Bases", len(cleaned_dna))
        metric_col2.metric("Codons", len(cleaned_dna) // 3)
        metric_col3.metric("GC content", f"{gc_content:.1f}%")
        chart_data = pd.DataFrame({"Base": list(counts), "Count": list(counts.values())})
        chart = alt.Chart(chart_data).mark_bar(color="#008f83").encode(x="Base", y="Count").properties(height=250)
        st.altair_chart(chart, use_container_width=True)
        if st.session_state.quality_passed:
            st.success("Quality check passed. The sequence is ready for transcription.")

with tab2:
    st.markdown('<div class="step-title">Step 2 · Transcription</div>', unsafe_allow_html=True)
    if not st.session_state.quality_passed:
        st.info("Complete the Quality Check first.")
    else:
        st.write("DNA thymine (T) is replaced by RNA uracil (U).")
        st.code(mrna_seq, language="text")
        if st.button("Approve transcription", key="approve_transcription"):
            st.session_state.transcription_passed = True
            st.rerun()
        if st.session_state.transcription_passed:
            st.success("Transcription approved. The mRNA is ready for translation.")

with tab3:
    st.markdown('<div class="step-title">Step 3 · Translation</div>', unsafe_allow_html=True)
    if not st.session_state.transcription_passed:
        st.info("Approve the transcription step first.")
    else:
        st.write("mRNA codons are converted into amino acids until the first stop codon.")
        st.code(protein_seq, language="text")
        st.metric("Amino acids", len(protein_seq))
        if not protein_seq:
            st.error("No protein was produced in this reading frame.")
        elif st.button("Approve protein", key="approve_translation"):
            st.session_state.translation_passed = True
            st.rerun()
        if st.session_state.translation_passed:
            st.success("Protein approved. It is ready for structure prediction.")

with tab4:
    st.markdown('<div class="step-title">Step 4 · Predict the 3D structure</div>', unsafe_allow_html=True)
    if not st.session_state.translation_passed:
        st.info("Approve the translation step first.")
    else:
        st.write("ESMFold will predict a structure from the validated amino-acid sequence.")
        st.code(protein_seq, language="text")
        if st.button("Generate 3D structure", key="generate_structure"):
            st.session_state.fold_error = None
            with st.spinner("Contacting ESMFold..."):
                pdb_data, error = get_3d_structure(protein_seq)
            st.session_state.pdb_data = pdb_data
            st.session_state.fold_error = error
            st.rerun()
        if st.session_state.fold_error:
            st.error(st.session_state.fold_error)
        if st.session_state.pdb_data:
            st.success("Structure generated.")
            view = py3Dmol.view(width=800, height=500)
            view.addModel(st.session_state.pdb_data, "pdb")
            view.setStyle({"cartoon": {"color": "spectrum"}})
            view.zoomTo()
            showmol(view, height=500, width=800)