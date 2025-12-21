import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Descriptors

# --- 1. SECURE CONFIGURATION ---
if "GEMINI_API_KEY" in st.secrets:
    genai.configure(api_key=st.secrets["GEMINI_API_KEY"])
    model = genai.GenerativeModel('gemini-1.5-flash') 
else:
    st.error("🔑 API Key Missing. Please add it to Streamlit Secrets.")
    st.stop()

# --- 2. THE STRATEGIC INVENTORY ---
# A list of the most problematic items for big-box retailers
product_inventory = {
    "Select a problematic item...": "",
    "Milk Jug (HDPE)": "CCCCCCCC",
    "Water Bottle (PET)": "CC1=CC=C(C=C1)C(=O)OCCO",
    "Aluminum Beverage Can": "[Al]",
    "Meat Wrap (PVC)": "C=CCl",
    "Deli Container (PP)": "CC(C)CC(C)C",
    "Frozen Food Bag (LDPE)": "CCCCCCCCCCCC",
    "Chip Bag (Multi-layer)": "CCCCCCCCCC.C=CC#N",
    "Styrofoam (EPS)": "C1=CC=C(C=C1)C=C",
    "Waxed Cardboard (Produce Box)": "CCCCCCCCCCCCCCCCCCCC",
    "Pharmacy Bottle (PC)": "CC(C)(C1=CC=C(OC(=O)OC2=CC=C(C(C)(C)C)C=C2)C=C1)C",
    "Coffee Cup Liner (PE)": "CCCCCCCC",
}

# --- 3. THE STRATEGIC DECISION ENGINE ---
def run_strategic_audit(item_name, smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol) if mol else 28.0
        toxic = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms()) if mol else False
        
        # PATH 1: RECYCLING (Resource Preservation)
        br = 92 if "PET" in item_name or "Aluminum" in item_name else max(10, 75 - (mw / 10))
        if toxic: br -= 30
        ar = min(99.1, br + 45) 
        
        # PATH 2: MINERALIZATION (Environmental Safety)
        bm = 5 if toxic else max(12, 50 - (mw / 20))
        am = 99.4 if not toxic else 94.0 
        
        # --- THE ARBITRATION LOGIC ---
        if toxic or "Multi-layer" in item_name or "Styrofoam" in item_name or "Waxed" in item_name:
            best_path = "Mineralization"
            priority = "Environmental Defense"
            reason = "This material is a 'Waste Trap.' Recycling requires toxic chemicals or excessive energy. Mineralization ensures zero-microplastic residue."
        elif ar > 90:
            best_path = "Mechanical Recycling"
            priority = "Resource Preservation"
            reason = "High-value circularity detected. Recovering this material saves more energy and carbon than returning it to the soil."
        else:
            best_path = "Mineralization"
            priority = "Environmental Defense"
            reason = "The efficiency gap for recycling is too high. Soil return is the safest environmental move."

        return round(br, 1), round(ar, 1), round(bm, 1), round(am, 1), best_path, priority, reason, mw
    except:
        return None

# --- 4. THE INTERFACE ---
st.set_page_config(page_title="VoraCycle Strategic Arbiter", layout="wide")
st.title("🔮 Wraith VoraCycle: Strategic Arbiter")
st.markdown("### *Forensic Search: Eliminating Environmental & Financial Debt*")

# DUAL SEARCH FUNCTIONALITY
col_search1, col_search2 = st.columns(2)
with col_search1:
    dropdown_query = st.selectbox("🧬 Preselected Problematic Items:", list(product_inventory.keys()))
with col_search2:
    manual_query = st.text_input("✍️ Manual Search (Type any item or material):")

# Logic to pick the query
query = manual_query if manual_query else dropdown_query

if query and query != "Select a problematic item...":
    smiles = product_inventory.get(query, "CCCCCCCC") # Defaults to generic polymer if not in list
    audit = run_strategic_audit(query, smiles)
    
    if audit:
        br, ar, bm, am, best_path, priority, reason, mw = audit
        
        st.divider()
        st.header(f"🏆 Best Strategic Path: {best_path}")
        st.markdown(f"**Primary Objective:** {priority}")
        st.info(f"**Decision Logic:** {reason}")
        
        # DUAL-PATH PERFORMANCE audit
        st.subheader("📊 Comparative Forensic Performance")
        p1, p2 = st.columns(2)
        with p1:
            st.metric("Mechanical Recycling (After)", f"{ar}%", delta=f"Before: {br}%")
            st.write("**Before Logic:** Low score due to **Chain Scission**. Heat-cycles degrade molecular integrity.")
            st.write("**After Logic:** High score via **Atomic Re-linkers** that heal bonds.")
        with p2:
            st.metric("Soil Mineralization (After)", f"{am}%", delta=f"Before: {bm}%")
            st.write("**Before Logic:** Low score due to **Biological Dead-Lock**. Microbes cannot 'see' the carbon.")
            st.write("**After Logic:** High score via **Metabolic Handles** that signal enzymes.")

        # --- FINAL FORENSIC CONCLUSION ---
        st.divider()
        st.header("📈 Strategic Forensic Deep-Dive")
        with st.spinner("Synthesizing Business & Environmental Impact..."):
            prompt = (
                f"Perform a deep forensic audit for {query}. The Arbiter chose {best_path} for {priority}. "
                f"1. Explain why the 'Before' state (Recycle: {br}%, Mineralize: {bm
