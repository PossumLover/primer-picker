import streamlit as st
from primers import create
from restriction_enzymes import RESTRICTION_ENZYMES
import pandas as pd
import json

# Try to import textcomplete for autocomplete; fall back if unavailable
try:
    from textcomplete import textcomplete, StrategyProps
    TEXTCOMPLETE_AVAILABLE = True
except ImportError:
    TEXTCOMPLETE_AVAILABLE = False
    st.warning("streamlit-textcomplete not installed. Install with: pip install streamlit-textcomplete")

# Pre-serialize your enzyme list for the JS search function
_ENZYME_LIST_JSON = json.dumps([
    {"name": name, "seq": seq}
    for name, seq in RESTRICTION_ENZYMES.items()
])

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
    Returns (results_text, forward_seq, reverse_seq)
    """
    if not target_sequence.strip():
        return "Error: Please provide a target DNA sequence.", "", ""
    
    # Clean up and validate the target
    seq = target_sequence.upper().replace(" ", "").replace("\n", "")
    if not set(seq).issubset(set("ATCG")):
        return "Error: Target sequence contains invalid characters. Only A, T, C, G are allowed.", "", ""
    
    # Build kwargs
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
    
    # Process any additions
    if add_fwd.strip():
        fwd_seq = process_enzyme_input(add_fwd.strip())
        if fwd_seq:
            kwargs['add_fwd'] = fwd_seq
    if add_rev.strip():
        rev_seq = process_enzyme_input(add_rev.strip())
        if rev_seq:
            kwargs['add_rev'] = rev_seq
    
    # Off-target check
    if parent_sequence.strip():
        parent = parent_sequence.upper().replace(" ", "").replace("\n", "")
        if not set(parent).issubset(set("ATCG")):
            return "Error: Parent sequence contains invalid characters. Only A, T, C, G are allowed.", "", ""
        kwargs['offtarget_check'] = parent
    
    try:
        fwd, rev = create(seq, **kwargs)
    except Exception as e:
        return f"Error designing primers: {e}", "", ""
    
    # Format results text
    results = f"""
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
    return results, fwd.seq, rev.seq

def process_enzyme_input(input_text: str) -> str:
    """
    Convert an enzyme name to its recognition sequence,
    or pass through a raw DNA sequence.
    """
    txt = input_text.strip()
    # Enzyme name?
    if txt in RESTRICTION_ENZYMES:
        return RESTRICTION_ENZYMES[txt]
    # Raw sequence?
    if set(txt.upper()).issubset(set("ATCG")):
        return txt.upper()
    # Otherwise, return as-is and let primers.create decide
    return txt

def setup_textcomplete_enzyme_input(label: str, key: str, placeholder: str = "") -> str:
    """
    Renders a single‐line text_input with optional
    streamlit-textcomplete autocomplete for enzyme names/sequences.
    Returns the final string (enzyme name or raw seq).
    """
    # 1) render the input box
    user_input = st.text_input(label, key=key, placeholder=placeholder)
    
    # 2) if autocomplete is available, hook it up
    if TEXTCOMPLETE_AVAILABLE:
        enzyme_strategy = StrategyProps(
            id=f"enzyme_{key}",
            # only trigger after 2+ word chars (letters, digits, underscore, hyphen)
            match=r"(\b[\w\-]{2,})$",
            # JS search filtering by name or seq prefix
            search=f"""
                async (term, callback) => {{
                  const enzymes = {_ENZYME_LIST_JSON};
                  const q = term.toLowerCase();
                  const results = enzymes
                    .filter(e =>
                      e.name.toLowerCase().startsWith(q) ||
                      e.seq.toLowerCase().startsWith(q)
                    )
                    .slice(0,10)
                    .map(e => `${{e.name}} (${{e.seq}})`);
                  callback(results);
                }}
            """,
            # on selection, replace the matched text with "EnzymeName "
            replace="value => value.split(' ')[0] + ' '",
            # simple HTML template for each dropdown item
            template="value => `<div style='display:flex;align-items:center'>🧬 ${value}</div>`"
        )
        textcomplete(
            area_label=label,
            strategies=[enzyme_strategy],
            max_count=10
        )
    
    # return whatever is in session_state for this key, stripped
    return st.session_state.get(key, "").strip()


# === Streamlit App ===

st.set_page_config(
    page_title="DNA Primer Designer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 DNA Primer Designer")
st.markdown("A simple site to design optimal PCR primers. Built by James Warner.")

# session state defaults
if 'target_seq_input' not in st.session_state:
    st.session_state.target_seq_input = ""
if 'parent_seq_input' not in st.session_state:
    st.session_state.parent_seq_input = ""

# layout
col1, col2 = st.columns([3,2], gap="large")
with col1:
    st.markdown("### Input Sequences")
    target_seq = st.text_area(
        "Target DNA Sequence (Required)",
        placeholder="A, T, C, G only",
        height=100,
        key="target_seq_input"
    )
    parent_seq = st.text_area(
        "Parent Sequence (Optional)",
        placeholder="For off-target checking",
        height=80,
        key="parent_seq_input"
    )
    st.markdown("### Primer Additions")
    if TEXTCOMPLETE_AVAILABLE:
        st.markdown("💡 Start typing enzyme names for autocomplete!")
    else:
        st.markdown("💡 Install streamlit-textcomplete for enzyme autocomplete.")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Forward Primer Addition**")
        add_fwd = setup_textcomplete_enzyme_input(
            "Forward Enzyme/Sequence",
            "fwd_enzyme_input",
            "e.g. BsaI or GGTCTC"
        )
        if add_fwd:
            st.caption(f"→ Will add: {process_enzyme_input(add_fwd)}")
    with c2:
        st.markdown("**Reverse Primer Addition**")
        add_rev = setup_textcomplete_enzyme_input(
            "Reverse Enzyme/Sequence",
            "rev_enzyme_input",
            "e.g. BpiI or GAAGAC"
        )
        if add_rev:
            st.caption(f"→ Will add: {process_enzyme_input(add_rev)}")

with col2:
    st.markdown("### Primer Parameters")
    st.markdown("**Optimal Values**")
    optimal_tm = st.slider("Optimal Tm (°C)", 50.0, 80.0, 62.0, 0.5)
    optimal_gc = st.slider("Optimal GC Content", 0.3, 0.7, 0.5, 0.05)
    optimal_len = st.slider("Optimal Length (bp)", 15, 35, 22, 1)
    st.markdown("**Penalty Weights**")
    penalty_tm = st.slider("Tm Penalty (per °C)", 0.1, 5.0, 1.0, 0.1)
    penalty_gc = st.slider("GC Penalty (%)", 0.1, 2.0, 0.2, 0.1)
    penalty_len = st.slider("Length Penalty (per bp)", 0.1, 2.0, 0.5, 0.1)
    penalty_tm_diff = st.slider("Tm Difference Penalty", 0.1, 5.0, 1.0, 0.1)
    penalty_dg = st.slider("Secondary Structure Penalty (ΔG)", 0.1, 10.0, 2.0, 0.1)

st.markdown("---")
if st.button("🔬 Design Primers", type="primary", use_container_width=True):
    if not target_seq.strip():
        st.error("Please enter a target DNA sequence.")
    else:
        with st.spinner("Designing primers..."):
            results_text, fwd_primer, rev_primer = design_primers(
                target_seq, parent_seq, add_fwd, add_rev,
                optimal_tm, optimal_gc, optimal_len,
                penalty_tm, penalty_gc, penalty_len,
                penalty_tm_diff, penalty_dg
            )
        st.markdown("### Results")
        c3, c4 = st.columns([2,1])
        with c3:
            st.text_area("Primer Design Results", value=results_text, height=400)
        with c4:
            if not results_text.startswith("Error"):
                st.text_input("Forward Primer (5' → 3')", value=fwd_primer)
                st.text_input("Reverse Primer (5' → 3')", value=rev_primer)
                df = pd.DataFrame({
                    'Primer': ['Forward','Reverse'],
                    'Sequence_5_to_3': [fwd_primer, rev_primer]
                })
                st.download_button(
                    "📥 Download Results as CSV",
                    df.to_csv(index=False),
                    file_name="primer_results.csv",
                    mime="text/csv"
                )

st.markdown("---")
st.markdown("### 📋 Example Sequences")
st.markdown("""
**Example 1:**  
Target: `AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA`  
Forward addition: `BsaI` → `GGTCTC`  
Reverse addition: `BpiI` → `GAAGAC`  

**Example 2 (no additions):**  
Target: `ATGAAACGCATTAGCACTGGGCCTAAGTACGAATTC`  

**Example 3 (off-target):**  
Target: `GCTAGCAATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAAGATCT`  
Parent: `ggaattacgtAATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAAggaccagttacagga`
""")

st.markdown("""
### 📖 Usage Instructions
1. Enter your **Target Sequence** (A/T/C/G only).  
2. (Optional) Enter a **Parent Sequence** for off-target checks.  
3. Under **Primer Additions**, start typing an enzyme name (e.g. “BsaI”) or a raw sequence (e.g. “GGTCTC”).  
   - If autocomplete is available, you’ll see suggestions like `BsaI (GGTCTC)`.  
   - Selecting a suggestion replaces your input with the enzyme name; we then convert it to the recognition site.  
4. Adjust **Optimal Values** and **Penalty Weights** as needed.  
5. Click **Design Primers** to generate your forward and reverse primers.  
6. Copy or download your results as CSV.
""")
