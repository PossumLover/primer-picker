import streamlit as st
from primers import create
from restriction_enzymes import RESTRICTION_ENZYMES
import pandas as pd
import json

# Try to import textcomplete, fall back to original implementation if not available
try:
    from textcomplete import textcomplete, StrategyProps
    TEXTCOMPLETE_AVAILABLE = True
except ImportError:
    TEXTCOMPLETE_AVAILABLE = False
    st.warning("streamlit-textcomplete not installed. Using fallback autocomplete. Install with: pip install streamlit-textcomplete")

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
        
        # Process forward primer addition
        if add_fwd.strip():
            fwd_sequence = process_enzyme_input(add_fwd.strip())
            if fwd_sequence:
                kwargs['add_fwd'] = fwd_sequence
        
        # Process reverse primer addition
        if add_rev.strip():
            rev_sequence = process_enzyme_input(add_rev.strip())
            if rev_sequence:
                kwargs['add_rev'] = rev_sequence
        
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

def process_enzyme_input(input_text):
    """
    Process enzyme input - convert enzyme name to sequence if needed.
    Returns the DNA sequence to be added to the primer.
    """
    input_text = input_text.strip().upper()
    
    # Check if input is an enzyme name
    if input_text in RESTRICTION_ENZYMES:
        return RESTRICTION_ENZYMES[input_text]
    
    # Check if input is already a DNA sequence
    valid_bases = set('ATCG')
    if all(base in valid_bases for base in input_text):
        return input_text
    
    # If not found and not a valid sequence, return as-is (let primers library handle)
    return input_text

def get_enzyme_suggestions(query):
    """Get enzyme suggestions based on name or sequence."""
    if not query:
        return []
    
    query = query.upper()
    suggestions = []
    
    # First try to match enzyme names (prioritize exact prefix matches)
    for enzyme, sequence in RESTRICTION_ENZYMES.items():
        if enzyme.upper().startswith(query):
            suggestions.append(f"{enzyme} ({sequence})")
    
    # Then try to match sequences
    for enzyme, sequence in RESTRICTION_ENZYMES.items():
        if sequence.startswith(query) and f"{enzyme} ({sequence})" not in suggestions:
            suggestions.append(f"{enzyme} ({sequence})")
    
    return suggestions[:10]  # Limit to 10 suggestions

def extract_enzyme_name_from_selection(selection):
    """Extract the enzyme name from the dropdown selection."""
    if not selection or '(' not in selection:
        return selection
    
    # Extract enzyme name from "EnzymeName (SEQUENCE)" format
    return selection.split(' (')[0]

def setup_textcomplete_enzyme_input(label, key, placeholder="e.g., BsaI or GGTCTC"):
    """Setup enzyme input with textcomplete autocomplete if available."""
    if TEXTCOMPLETE_AVAILABLE:
        # Create text input instead of text area for better autocomplete experience
        enzyme_input = st.text_input(
            label=label,
            placeholder=placeholder,
            key=key,
            max_chars=50
        )
        
        # Create the search function that returns enzyme suggestions
        def create_search_function():
            return f"""
            async (term, callback) => {{
                const enzymes = {json.dumps([
                    {"name": enzyme, "sequence": sequence, "display": f"{enzyme} ({sequence})"}
                    for enzyme, sequence in RESTRICTION_ENZYMES.items()
                ])};
                
                const query = term.toLowerCase();
                const results = enzymes
                    .filter(enzyme => 
                        enzyme.name.toLowerCase().startsWith(query) || 
                        enzyme.sequence.toLowerCase().startsWith(query)
                    )
                    .slice(0, 10)
                    .map(enzyme => enzyme.display);
                
                callback(results);
            }}
            """
        
        # Create autocomplete strategy
        enzyme_strategy = StrategyProps(
            id=f"enzyme_{key}",
            match=r"(\w*)$",  # Match word characters at end of input
            search=create_search_function(),
            replace="(value) => value.split(' ')[0]",  # Extract just the enzyme name
            template="(value) => `🧬 ${value}`",
        )
        
        # Initialize textcomplete with proper configuration
        textcomplete(
            area_label=label,
            strategies=[enzyme_strategy],
            max_count=10,
            dropdown_className="enzyme-autocomplete-dropdown"
        )
        
        return enzyme_input.strip() if enzyme_input else ""
    else:
        # Fallback to original implementation with improved UX
        enzyme_input = st.text_input(
            label,
            placeholder=placeholder,
            key=key
        )
        
        # Show live suggestions if there's input
        if enzyme_input and len(enzyme_input) >= 2:  # Only show suggestions after 2+ characters
            suggestions = get_enzyme_suggestions(enzyme_input)
            if suggestions:
                # Create a more user-friendly selection
                selected = st.selectbox(
                    "💡 Suggestions (or continue typing):",
                    options=["Keep typing..."] + suggestions,
                    key=f"{key}_select",
                    help="Select a suggestion or continue typing your own sequence"
                )
                
                if selected == "Keep typing...":
                    return enzyme_input
                else:
                    # Extract enzyme name and update the input
                    enzyme_name = extract_enzyme_name_from_selection(selected)
                    # Use st.session_state to update the input value
                    st.session_state[key] = enzyme_name
                    return enzyme_name
            else:
                # Show helpful message when no matches found
                if all(c in 'ATCG' for c in enzyme_input.upper()):
                    st.info("✅ Custom DNA sequence detected")
                else:
                    st.info("🔍 No matching enzymes found. Try typing an enzyme name or DNA sequence.")
                return enzyme_input
        
        return enzyme_input

# Streamlit app
st.set_page_config(
    page_title="DNA Primer Designer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 DNA Primer Designer")
st.markdown("A simple website to help design optimal PCR primers for your DNA sequences. Built by James Warner.")

# Initialize session state variables if they don't exist
if 'target_seq_input' not in st.session_state:
    st.session_state.target_seq_input = ""
if 'parent_seq_input' not in st.session_state:
    st.session_state.parent_seq_input = ""

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
    
    if TEXTCOMPLETE_AVAILABLE:
        st.markdown("💡 **Pro tip**: Start typing enzyme names for autocomplete suggestions!")
    else:
        st.markdown("💡 **Tip**: Type enzyme names (e.g., 'BsaI') and select from suggestions, or enter DNA sequences directly.")
    
    # Forward and reverse primer additions with autocomplete
    col1a, col1b = st.columns(2)
    
    with col1a:
        st.markdown("**Forward Primer Addition**")
        add_fwd = setup_textcomplete_enzyme_input(
            "Forward Enzyme/Sequence",
            "fwd_enzyme_input",
            "e.g., BsaI or GGTCTC"
        )
        
        # Show what sequence will be added
        if add_fwd:
            processed_seq = process_enzyme_input(add_fwd)
            if processed_seq != add_fwd.upper():
                st.caption(f"Will add: {processed_seq}")
    
    with col1b:
        st.markdown("**Reverse Primer Addition**")
        add_rev = setup_textcomplete_enzyme_input(
            "Reverse Enzyme/Sequence", 
            "rev_enzyme_input",
            "e.g., BpiI or GAAGAC"
        )
        
        # Show what sequence will be added
        if add_rev:
            processed_seq = process_enzyme_input(add_rev)
            if processed_seq != add_rev.upper():
                st.caption(f"Will add: {processed_seq}")

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

# Examples section for reference
st.markdown("---")
st.markdown("### 📋 Example Sequences")
st.markdown("""
**Example 1: Basic sequence with BsaI/BpiI sites**
- Target: `AATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAA`
- Forward addition: `GGTCTC` (BsaI) or just type "BsaI"
- Reverse addition: `GAAGAC` (BpiI) or just type "BpiI"

**Example 2: Simple sequence without additions**
- Target: `ATGAAACGCATTAGCACTGGGCCTAAGTACGAATTC`

**Example 3: With parent sequence for off-target checking**
- Target: `GCTAGCAATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAAGATCT`
- Parent: `ggaattacgtAATGAGACAATAGCACACACAGCTAGGTCAGCATACGAAAggaccagttacagga`
""")

# Instructions
st.markdown("""
### 📖 Usage Instructions:

1. **Target Sequence**: Enter your DNA sequence of interest (required)
2. **Parent Sequence**: Optionally provide a larger sequence context to check for off-target binding
3. **Primer Additions**: Add restriction enzyme sites or other sequences to your primers
   - **With autocomplete**: Start typing enzyme names (e.g., "Bsa") and select from suggestions
   - **Fallback mode**: Type enzyme names (e.g., "BsaI") or sequences (e.g., "GGTCTC") directly
   - The system will automatically convert enzyme names to their recognition sequences
4. **Optimal Parameters**: Set target values for primer characteristics
5. **Penalty Weights**: Adjust how strongly deviations from optimal are penalized
6. **Click "Design Primers"** to generate optimal forward and reverse primers

**How Enzyme Input Works:**
- Type enzyme names like "BsaI", "BpiI", "EcoRI" - they'll be converted to recognition sequences
- Or enter DNA sequences directly like "GGTCTC", "GAAGAC"
- Autocomplete suggestions show both enzyme name and sequence for reference

**Penalty Calculation**: Each primer's penalty = |Tm deviation| × Tm penalty + |GC deviation| × GC penalty + |length deviation| × length penalty + |primer Tm difference| × Tm diff penalty + |ΔG| × ΔG penalty

The algorithm selects primers with the lowest total penalty scores.
""")
