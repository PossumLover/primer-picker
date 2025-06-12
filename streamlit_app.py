import streamlit as st
from primers import create
import pandas as pd

# Restriction enzyme database
RESTRICTION_ENZYMES = {
    'AclI': 'AACGTT',
    'HindIII': 'AAGCTT',
    'SspI-HF': 'AATATT',
    'MluCI': 'AATT',
    'PciI': 'ACATGT',
    'AgeI-HF': 'ACCGGT',
    'BfuAI': 'ACCTGC',
    'BspMI': 'ACCTGC',
    'SexAI': 'ACCWGGT',
    'MluI-HF': 'ACGCGT',
    'BceAI': 'ACGGC',
    'HpyCH4IV': 'ACGT',
    'HpyCH4III': 'ACNGT',
    'BaeI': 'ACNNNNGTAYC',
    'BsaXI': 'ACNNNNNCTCC',
    'AflIII': 'ACRYGT',
    'SpeI-HF': 'ACTAGT',
    'BsrI': 'ACTGG',
    'BmrI': 'ACTGGG',
    'BglII': 'AGATCT',
    'AfeI': 'AGCGCT',
    'AluI': 'AGCT',
    'StuI': 'AGGCCT',
    'ScaI-HF': 'AGTACT',
    'ClaI': 'ATCGAT',
    'BspDI': 'ATCGAT',
    'PI-SceI': 'ATCTATGTCGGGTGCGGAGAAAGAGGTAAT',
    'NsiI-HF': 'ATGCAT',
    'NsiI': 'ATGCAT',
    'AseI': 'ATTAAT',
    'SwaI': 'ATTTAAAT',
    'CspCI': 'CAANNNNNGTGG',
    'MfeI-HF': 'CAATTG',
    'PaqCI': 'CACCTGC',
    'BssSI-v2': 'CACGAG',
    'Nb.BssSI': 'CACGAG',
    'BmgBI': 'CACGTC',
    'PmlI': 'CACGTG',
    'DraIII-HF': 'CACNNNGTG',
    'AleI-v2': 'CACNNNNGTG',
    'EcoP15I': 'CAGCAG',
    'PvuII-HF': 'CAGCTG',
    'PvuII': 'CAGCTG',
    'AlwNI': 'CAGNNNCTG',
    'BtsIMutI': 'CAGTG',
    'NdeI': 'CATATG',
    'NlaIII': 'CATG',
    'FatI': 'CATG',
    'MslI': 'CAYNNNNRTG',
    'XcmI': 'CCANNNNNNNNNTGG',
    'BstXI': 'CCANNNNNNTGG',
    'PflMI': 'CCANNNNNTGG',
    'BccI': 'CCATC',
    'NcoI-HF': 'CCATGG',
    'NcoI': 'CCATGG',
    'BseYI': 'CCCAGC',
    'FauI': 'CCCGC',
    'TspMI': 'CCCGGG',
    'XmaI': 'CCCGGG',
    'SmaI': 'CCCGGG',
    'Nt.CviPII': 'CCD',
    'AciI': 'CCGC',
    'SacII': 'CCGCGG',
    'BsrBI': 'CCGCTC',
    'HpaII': 'CCGG',
    'MspI': 'CCGG',
    'ScrFI': 'CCNGG',
    'StyD4I': 'CCNGG',
    'BsaJI': 'CCNNGG',
    'BslI': 'CCNNNNNNNGG',
    'BtgI': 'CCRYGG',
    'NciI': 'CCSGG',
    'AvrII': 'CCTAGG',
    'MnlI': 'CCTC',
    'BbvCI': 'CCTCAGC',
    'Nt.BbvCI': 'CCTCAGC',
    'Nb.BbvCI': 'CCTCAGC',
    'SbfI-HF': 'CCTGCAGG',
    'Bpu10I': 'CCTNAGC',
    'Bsu36I': 'CCTNAGG',
    'EcoNI': 'CCTNNNNNAGG',
    'HpyAV': 'CCTTC',
    'BstNI': 'CCWGG',
    'PspGI': 'CCWGG',
    'StyI-HF': 'CCWWGG',
    'BcgI': 'CGANNNNNNTGC',
    'PvuI-HF': 'CGATCG',
    'BstUI': 'CGCG',
    'EagI-HF': 'CGGCCG',
    'RsrII': 'CGGWCCG',
    'BsiEI': 'CGRYCG',
    'BsiWI-HF': 'CGTACG',
    'BsiWI': 'CGTACG',
    'BsmBI-v2': 'CGTCTC',
    'Esp3I': 'CGTCTC',
    'Hpy99I': 'CGWCG',
    'MspA1I': 'CMGCKG',
    'MspJI': 'CNNR',
    'SgrAI': 'CRCCGGYG',
    'BfaI': 'CTAG',
    'BspCNI': 'CTCAG',
    'XhoI': 'CTCGAG',
    'PaeR7I': 'CTCGAG',
    'EarI': 'CTCTTC',
    'AcuI': 'CTGAAG',
    'PstI-HF': 'CTGCAG',
    'PstI': 'CTGCAG',
    'BpmI': 'CTGGAG',
    'DdeI': 'CTNAG',
    'SfcI': 'CTRYAG',
    'AflII': 'CTTAAG',
    'BpuEI': 'CTTGAG',
    'SmlI': 'CTYRAG',
    'BsoBI': 'CYCGRG',
    'AvaI': 'CYCGRG',
    'MboII': 'GAAGA',
    'BbsI': 'GAAGAC',
    'BbsI-HF': 'GAAGAC',
    'XmnI': 'GAANNNNTTC',
    'BsmI': 'GAATGC',
    'Nb.BsmI': 'GAATGC',
    'EcoRI': 'GAATTC',
    'EcoRI-HF': 'GAATTC',
    'HgaI': 'GACGC',
    'AatII': 'GACGTC',
    'ZraI': 'GACGTC',
    'PflFI': 'GACNNNGTC',
    'Tth111I': 'GACNNNGTC',
    'PshAI': 'GACNNNNGTC',
    'AhdI': 'GACNNNNNGTC',
    'DrdI': 'GACNNNNNNGTC',
    'Eco53kI': 'GAGCTC',
    'SacI-HF': 'GAGCTC',
    'BseRI': 'GAGGAG',
    'MlyI': 'GAGTC',
    'PleI': 'GAGTC',
    'Nt.BstNBI': 'GAGTC',
    'HinfI': 'GANTC',
    'EcoRV': 'GATATC',
    'EcoRV-HF': 'GATATC',
    'DpnI': 'GATC',
    'MboI': 'GATC',
    'Sau3AI': 'GATC',
    'DpnII': 'GATC',
    'BsaBI': 'GATNNNNNATC',
    'TfiI': 'GAWTC',
    'BsrDI': 'GCAATG',
    'Nb.BsrDI': 'GCAATG',
    'BbvI': 'GCAGC',
    'BtsI-v2': 'GCAGTG',
    'Nb.BtsI': 'GCAGTG',
    'BstAPI': 'GCANNNNNTGC',
    'SfaNI': 'GCATC',
    'SphI-HF': 'GCATGC',
    'SphI': 'GCATGC',
    'SrfI': 'GCCCGGGC',
    'NmeAIII': 'GCCGAG',
    'NgoMIV': 'GCCGGC',
    'NaeI': 'GCCGGC',
    'BglI': 'GCCNNNNNGCC',
    'AsiSI': 'GCGATCGC',
    'BtgZI': 'GCGATG',
    'HinP1I': 'GCGC',
    'HhaI': 'GCGC',
    'BssHII': 'GCGCGC',
    'NotI-HF': 'GCGGCCGC',
    'NotI': 'GCGGCCGC',
    'Fnu4HI': 'GCNGC',
    'Cac8I': 'GCNNGC',
    'MwoI': 'GCNNNNNNNCC',
    'BmtI-HF': 'GCTAGC',
    'NheI-HF': 'GCTAGC',
    'BspQI-HF': 'GCTCTTC',
    'BspQI': 'GCTCTTC',
    'SapI': 'GCTCTTC',
    'Nt.BspQI': 'GCTCTTC',
    'BlpI': 'GCTNAGC',
    'ApeKI': 'GCWGC',
    'TseI': 'GCWGC',
    'Bsp1286I': 'GDGCHC',
    'AlwI': 'GGATC',
    'BamHI-HF': 'GGATCC',
    'BamHI': 'GGATCC',
    'Nt.AlwI': 'GGATC',
    'FokI': 'GGATG',
    'BtsCI': 'GGATG',
    'HaeIII': 'GGCC',
    'FseI': 'GGCCGGCC',
    'SfiI': 'GGCCNNNNNGGCC',
    'NarI': 'GGCGCC',
    'SfoI': 'GGCGCC',
    'PluTI': 'GGCGCC',
    'KasI': 'GGCGCC',
    'AscI': 'GGCGCGCC',
    'EciI': 'GGCGGA',
    'BsmFI': 'GGGAC',
    'ApaI': 'GGGCCC',
    'PspOMI': 'GGGCCC',
    'Sau96I': 'GGNCC',
    'NlaIV': 'GGNNCC',
    'KpnI-HF': 'GGTACC',
    'Acc65I': 'GGTACC',
    'BsaI-HF': 'GGTCTC',
    'HphI': 'GGTGA',
    'BstEII-HF': 'GGTNACC',
    'AvaII': 'GGWCC',
    'BanI': 'GGYRCC',
    'BaeGI': 'GKGCMC',
    'BsaHI': 'GRCGYC',
    'BanII': 'GRGCYC',
    'CviQI': 'GTAC',
    'RsaI': 'GTAC',
    'BstZ17I-HF': 'GTATAC',
    'BciVI': 'GTATCC',
    'SalI': 'GTCGAC',
    'SalI-HF': 'GTCGAC',
    'BsmAI': 'GTCTC',
    'BcoDI': 'GTCTC',
    'Nt.BsmAI': 'GTCTC',
    'ApaLI': 'GTGCAC',
    'BsgI': 'GTGCAG',
    'AccI': 'GTMKAC',
    'Hpy166II': 'GTNNAC',
    'Tsp45I': 'GTSAC',
    'HpaI': 'GTTAAC',
    'PmeI': 'GTTTAAAC',
    'HincII': 'GTYRAC',
    'BsiHKAI': 'GWGCWC',
    'TspRI': 'NNCASTGNN',
    'ApoI-HF': 'RAATTY',
    'NspI': 'RCATGY',
    'BsrFI-v2': 'RCCGGY',
    'BstYI': 'RGATCY',
    'HaeII': 'RGCGCY',
    'CviKI-1': 'RGCY',
    'EcoO109I': 'RGGNCCY',
    'PpuMI': 'RGGWCCY',
    'I-CeuI': 'TAACTATAACGGTCCTAAGGTAGCGAA',
    'SnaBI': 'TACGTA',
    'I-SceI': 'TAGGGATAACAGGGTAAT',
    'BspHI': 'TCATGA',
    'BspEI': 'TCCGGA',
    'MmeI': 'TCCRAC',
    'TaqI-v2': 'TCGA',
    'NruI-HF': 'TCGCGA',
    'Hpy188I': 'TCNGA',
    'Hpy188III': 'TCNNGA',
    'XbaI': 'TCTAGA',
    'BclI': 'TGATCA',
    'BclI-HF': 'TGATCA',
    'HpyCH4V': 'TGCA',
    'FspI': 'TGCGCA',
    'PI-PspI': 'TGGCAAACAGCTATTATGGGTATTATGGGT',
    'MscI': 'TGGCCA',
    'BsrGI-HF': 'TGTACA',
    'MseI': 'TTAA',
    'PacI': 'TTAATTAA',
    'PsiI-v2': 'TTATAA',
    'BstBI': 'TTCGAA',
    'DraI': 'TTTAAA',
    'PspXI': 'VCTCGAGB',
    'BsaWI': 'WCCGGW',
    'BsaAI': 'YACGTR',
    'EaeI': 'YGGCCR'
}

def design_primers(
    target_sequence,
    parent_sequence="",
    add_fwd="",
    add_rev="",
    optimal_tm=62.0,
    optimal_gc=0.5,
    optimal_len=22,
    penalty_tm=1.0,
    penalty_gc=0.2,
    penalty_len=0.5,
    penalty_tm_diff=1.0,
    penalty_dg=2.0
):
    """
    Design primers for a given DNA sequence using the primers library.
    """
    if not target_sequence.strip():
        return "Error: Please provide a target DNA sequence.", "", ""
    
    # Clean up the sequence (remove whitespace and convert to uppercase)
    target_sequence = target_sequence.strip().upper().replace(" ", "").replace("\n", "")
    
    # Validate DNA sequence
    valid_bases = set('ATCG')
    if not all(base in valid_bases for base in target_sequence):
        return "Error: Target sequence contains invalid characters. Only A, T, C, G are allowed.", "", ""
    
    try:
        # Prepare arguments for primers.create()
        kwargs = {
            'optimal_tm': optimal_tm,
            'optimal_gc': optimal_gc,
            'optimal_len': int(optimal_len),
            'penalty_tm': penalty_tm,
            'penalty_gc': penalty_gc,
            'penalty_len': penalty_len,
            'penalty_tm_diff': penalty_tm_diff,
            'penalty_dg': penalty_dg
        }
        
        # Add forward primer additions if specified
        if add_fwd.strip():
            kwargs['add_fwd'] = add_fwd.strip().upper()
        
        # Add reverse primer additions if specified
        if add_rev.strip():
            kwargs['add_rev'] = add_rev.strip().upper()
        
        # Add parent sequence for off-target checking if specified
        if parent_sequence.strip():
            parent_sequence = parent_sequence.strip().upper().replace(" ", "").replace("\n", "")
            # Validate parent sequence
            if not all(base in valid_bases for base in parent_sequence):
                return "Error: Parent sequence contains invalid characters. Only A, T, C, G are allowed.", "", ""
            kwargs['offtarget_check'] = parent_sequence
        
        # Create primers
        fwd, rev = create(target_sequence, **kwargs)
        
        # Format results
        results_text = f"""
PRIMER DESIGN RESULTS
====================

FORWARD PRIMER:
Sequence (5' to 3'): {fwd.seq}
Length: {len(fwd.seq)} bp
Melting Temperature (Tm): {fwd.tm:.1f}°C
Total Tm (with additions): {fwd.tm_total:.1f}°C
GC Content: {fwd.gc:.1%}
Free Energy (ΔG): {fwd.dg:.2f} kcal/mol

REVERSE PRIMER:
Sequence (5' to 3'): {rev.seq}
Length: {len(rev.seq)} bp
Melting Temperature (Tm): {rev.tm:.1f}°C
Total Tm (with additions): {rev.tm_total:.1f}°C
GC Content: {rev.gc:.1%}
Free Energy (ΔG): {rev.dg:.2f} kcal/mol

PRIMER PAIR ANALYSIS:
Tm Difference: {abs(fwd.tm_total - rev.tm_total):.1f}°C
"""
        
        # Create CSV data for download
        primer_data = {
            'Primer': ['Forward', 'Reverse'],
            'Sequence_5_to_3': [fwd.seq, rev.seq],
            'Length_bp': [len(fwd.seq), len(rev.seq)],
            'Tm_C': [fwd.tm, rev.tm],
            'Total_Tm_C': [fwd.tm_total, rev.tm_total],
            'GC_Content': [fwd.gc, rev.gc],
            'Delta_G_kcal_mol': [fwd.dg, rev.dg]
        }
        
        df = pd.DataFrame(primer_data)
        csv_output = df.to_csv(index=False)
        
        return results_text, fwd.seq, rev.seq
        
    except Exception as e:
        return f"Error designing primers: {str(e)}", "", ""

def get_enzyme_suggestions(query):
    """Get enzyme suggestions based on name or sequence."""
    if not query:
        return []
    
    query = query.upper()
    suggestions = []
    
    # First try to match enzyme names
    for enzyme, sequence in RESTRICTION_ENZYMES.items():
        if enzyme.upper().startswith(query):
            suggestions.append(f"{enzyme} ({sequence})")
    
    # Then try to match sequences
    for enzyme, sequence in RESTRICTION_ENZYMES.items():
        if sequence.startswith(query) and f"{enzyme} ({sequence})" not in suggestions:
            suggestions.append(f"{enzyme} ({sequence})")
    
    return suggestions[:10]  # Limit to 10 suggestions

def extract_sequence_from_selection(selection):
    """Extract the recognition sequence from the dropdown selection."""
    if not selection or '(' not in selection:
        return selection
    
    # Extract sequence from "EnzymeName (SEQUENCE)" format
    start = selection.find('(') + 1
    end = selection.find(')')
    if start > 0 and end > start:
        return selection[start:end]
    return selection

# Streamlit app
st.set_page_config(
    page_title="DNA Primer Designer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 DNA Primer Designer")
st.markdown("Design optimal PCR primers for your DNA sequences using the `primers` Python library.")

# Create main layout with proper spacing
col1, col2 = st.columns([3, 2], gap="large")

with col1:
    st.markdown("### Input Sequences")
    
    target_seq = st.text_area(
        "Target DNA Sequence (Required)",
        placeholder="Enter your target DNA sequence (A, T, C, G only)",
        height=100,
        key="target_seq_input"
    )
    
    parent_seq = st.text_area(
        "Parent Sequence (Optional - for off-target checking)",
        placeholder="Enter parent sequence to check for off-target binding",
        height=80,
        key="parent_seq_input"
    )
    
    st.markdown("### Primer Additions")
    st.markdown("*Add restriction enzyme recognition sites or custom sequences to your primers*")
    
    # Forward primer addition with smart enzyme selection
    col1a, col1b = st.columns(2)
    
    with col1a:
        st.markdown("**Forward Primer Addition**")
        fwd_enzyme_input = st.text_input(
            "Type enzyme name or sequence",
            placeholder="e.g., BsaI or GGTCTC",
            key="fwd_enzyme_search",
            label_visibility="collapsed"
        )
        
        # Get suggestions for forward primer
        fwd_suggestions = get_enzyme_suggestions(fwd_enzyme_input)
        if fwd_suggestions:
            fwd_selected = st.selectbox(
                "Select enzyme (Forward)",
                options=[""] + fwd_suggestions,
                key="fwd_enzyme_select",
                label_visibility="collapsed"
            )
            add_fwd = extract_sequence_from_selection(fwd_selected) if fwd_selected else fwd_enzyme_input
        else:
            add_fwd = fwd_enzyme_input
    
    with col1b:
        st.markdown("**Reverse Primer Addition**")
        rev_enzyme_input = st.text_input(
            "Type enzyme name or sequence",
            placeholder="e.g., BpiI or GAAGAC",
            key="rev_enzyme_search",
            label_visibility="collapsed"
        )
        
        # Get suggestions for reverse primer
        rev_suggestions = get_enzyme_suggestions(rev_enzyme_input)
        if rev_suggestions:
            rev_selected = st.selectbox(
                "Select enzyme (Reverse)",
                options=[""] + rev_suggestions,
                key="rev_enzyme_select",
                label_visibility="collapsed"
            )
            add_rev = extract_sequence_from_selection(rev_selected) if rev_selected else rev_enzyme_input
        else:
            add_rev = rev_enzyme_input

with col2:
    st.markdown("### Primer Parameters")
    
    st.markdown("**Optimal Values:**")
    optimal_tm = st.slider(
        "Optimal Tm (°C)",
        min_value=50.0,
        max_value=80.0,
        value=62.0,
        step=0.5
    )
    
    optimal_gc = st.slider(
        "Optimal GC Content",
        min_value=0.3,
        max_value=0.7,
        value=0.5,
        step=0.05
    )
    
    optimal_len = st.slider(
        "Optimal Length (bp)",
        min_value=15,
        max_value=35,
        value=22,
        step=1
    )
    
    st.markdown("**Penalty Weights:**")
    st.markdown("*Higher values = stronger penalties for deviations from optimal*")
    
    penalty_tm = st.slider(
        "Tm Penalty (per °C deviation)",
        min_value=0.1,
        max_value=5.0,
        value=1.0,
        step=0.1,
        help="Penalizes melting temperature deviation from optimal"
    )
    
    penalty_gc = st.slider(
        "GC Content Penalty (per % deviation)",
        min_value=0.1,
        max_value=2.0,
        value=0.2,
        step=0.1,
        help="Penalizes GC content deviation from optimal"
    )
    
    penalty_len = st.slider(
        "Length Penalty (per bp deviation)",
        min_value=0.1,
        max_value=2.0,
        value=0.5,
        step=0.1,
        help="Penalizes primer length deviation from optimal"
    )
    
    penalty_tm_diff = st.slider(
        "Tm Difference Penalty (per °C diff between primers)",
        min_value=0.1,
        max_value=5.0,
        value=1.0,
        step=0.1,
        help="Penalizes melting temperature difference between forward and reverse primers"
    )
    
    penalty_dg = st.slider(
        "Secondary Structure Penalty (per kcal/mol)",
        min_value=0.1,
        max_value=10.0,
        value=2.0,
        step=0.1,
        help="Penalizes free energy in secondary structures (higher |ΔG| = more stable structures)"
    )

# Design button with better spacing
st.markdown("---")
if st.button("🔬 Design Primers", type="primary", use_container_width=True):
    if target_seq.strip():
        with st.spinner("Designing primers..."):
            results_text, fwd_primer, rev_primer = design_primers(
                target_seq, parent_seq, add_fwd, add_rev,
                optimal_tm, optimal_gc, optimal_len,
                penalty_tm, penalty_gc, penalty_len, penalty_tm_diff, penalty_dg
            )
        
        st.markdown("### Results")
        
        # Display results
        col3, col4 = st.columns([2, 1])
        
        with col3:
            st.text_area(
                "Primer Design Results",
                value=results_text,
                height=400
            )
        
        with col4:
            if not results_text.startswith("Error"):
                st.text_input(
                    "Forward Primer (5' → 3')",
                    value=fwd_primer
                )
                st.text_input(
                    "Reverse Primer (5' → 3')",
                    value=rev_primer
                )
                
                # Create downloadable CSV
                primer_data = {
                    'Primer': ['Forward', 'Reverse'],
                    'Sequence_5_to_3': [fwd_primer, rev_primer]
                }
                df = pd.DataFrame(primer_data)
                csv = df.to_csv(index=False)
                
                st.download_button(
                    label="📥 Download Results as CSV",
                    data=csv,
                    file_name="primer_results.csv",
                    mime="text/csv"
                )
    else:
        st.error("Please enter a target DNA sequence.")

# Examples section with better formatting
st.markdown("---")
st.markdown("### 📋 Example Sequences")

examples = [
    {
        "name": "Example 1: Basic sequence with BsaI/BpiI sites",
        "target": "AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA",
        "parent": "",
        "fwd_add": "GGTCTC",
        "rev_add": "GAAGAC"
    },
    {
        "name": "Example 2: Simple sequence without additions",
        "target": "ATGAAACGCATTAGCACTGGGCCTAAGTACGAATTC",
        "parent": "",
        "fwd_add": "",
        "rev_add": ""
    },
    {
        "name": "Example 3: With parent sequence for off-target checking",
        "target": "GCTAGCAATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAAGATCT",
        "parent": "ggaattacgtAATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAAggaccagttacagga",
        "fwd_add": "",
        "rev_add": ""
    }
]

selected_example = st.selectbox(
    "Choose an example to load:",
    options=[""] + [ex["name"] for ex in examples],
    key="example_selector"
)

if selected_example:
    example = next(ex for ex in examples if ex["name"] == selected_example)
    if st.button("Load Example", key="load_example"):
        st.session_state.target_seq = example["target"]
        st.session_state.parent_seq = example["parent"]
        st.session_state.add_fwd = example["fwd_add"]
        st.session_state.add_rev = example["rev_add"]
        st.rerun()

# Instructions
st.markdown("""
### 📖 Usage Instructions:

1. **Target Sequence**: Enter your DNA sequence of interest (required)
2. **Parent Sequence**: Optionally provide a larger sequence context to check for off-target binding
3. **Primer Additions**: Add restriction enzyme sites or other sequences to your primers
4. **Optimal Parameters**: Set target values for primer characteristics
5. **Penalty Weights**: Adjust how strongly deviations from optimal are penalized
6. **Click "Design Primers"** to generate optimal forward and reverse primers

**Penalty Calculation**: Each primer's penalty = |Tm deviation| × Tm penalty + |GC deviation| × GC penalty + |length deviation| × length penalty + |primer Tm difference| × Tm diff penalty + |ΔG| × ΔG penalty

The algorithm selects primers with the lowest total penalty scores.
""")
