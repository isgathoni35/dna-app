import streamlit as st
import pandas as pd
import altair as alt

# 1. Page Configuration
st.set_page_config(page_title="DNA Nucleotide Count", page_icon="🧬")

# 2. App Title and Description
st.write("""
# 🧬 DNA Nucleotide Count Web App
This app counts the composition of nucleotides (A, T, G, C) in your query DNA sequence.
***
""")

# 3. Sidebar Input
st.sidebar.header('Enter DNA Sequence')
sequence_input = ">DNA Query\nGAACACGTGGAGGCAAACAGGAAGGTGAAGAAGAACTTATCCTATCAGGACGGAAGGTCCTGTGCTCGGG\nATCTTCCAGACGTCGCGACTCTAAATTGCCCCCTCTGAGGTCAAGGAACACAAGATGGTTTTGGAAATGC\nTGAACCCGATACATTATAACATCACCAGCATCGTGCCTGAAGCCATGCCTGCTGCCACCATGCCAGTCCT"

sequence = st.sidebar.text_area("Sequence Input", sequence_input, height=250)
sequence = sequence.splitlines()
sequence = sequence[1:] # Skips the first line (the name of the sequence)
sequence = ''.join(sequence) # Concatenates the list to string

st.write("""
***
""")

# 4. Display the Input
st.header('INPUT (DNA Query)')
st.code(sequence)

# 5. Calculate Nucleotides
st.header('OUTPUT (DNA Composition)')

# Custom function to count
def DNA_nucleotide_count(seq):
  d = dict([
            ('A', seq.count('A')),
            ('T', seq.count('T')),
            ('G', seq.count('G')),
            ('C', seq.count('C'))
            ])
  return d

X = DNA_nucleotide_count(sequence)

# 6. Display Counts
st.subheader('1. Print Dictionary')
st.write(X)

st.subheader('2. Print Text')
st.write('There are  ' + str(X['A']) + ' adenine (A)')
st.write('There are  ' + str(X['T']) + ' thymine (T)')
st.write('There are  ' + str(X['G']) + ' guanine (G)')
st.write('There are  ' + str(X['C']) + ' cytosine (C)')

# 7. Display DataFrame
st.subheader('3. Display DataFrame')
df = pd.DataFrame.from_dict(X, orient='index')
df = df.rename({0: 'count'}, axis='columns')
df.reset_index(inplace=True)
df = df.rename(columns = {'index':'nucleotide'})
st.write(df)

# 8. Display Bar Chart
st.subheader('4. Display Bar chart')
p = alt.Chart(df).mark_bar().encode(
    x='nucleotide',
    y='count'
)
p = p.properties(
    width=alt.Step(80)  # controls width of bar.
)
st.write(p)