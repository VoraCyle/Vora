import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Descriptors

# --- 1. SECURE CONFIGURATION ---
# The logic gate: If the key is found, configure the AI. If not, stop.
if "GEMINI_API_KEY" in st.secrets:
    genai.configure(api_key=st.secrets "AIzaSyBEmthiZ4aKWUONVhvM4XAU9dCofarQ6EQ")
    model = genai.GenerativeModel('gemini-1.5-flash') 
else:

# --- 2. THE STRATEGIC INVENTORY ---
# This dictionary must be at the ROOT level (zero spaces at the start)
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
}

# --- 3. THE STRATEGIC DECISION ENGINE ---
def run_strategic_audit(item_name, smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol) if mol else 28.0
        toxic = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms()) if mol else False
        
        # Before scores represent the baseline liability
        br = 92 if "PET" in item_name or "Aluminum" in item_name else max(10, 75 - (mw / 10))
        if toxic: br -= 30
        ar = min(99.1, br + 45) # After score via VoraCycle Surgery
        
        bm = 5 if toxic else max(12, 50 - (mw / 20))
        am = 99.4 if not toxic else 94.0 
        
        # Determining the Apex Move
        if toxic or "Multi-layer" in item_name or "Styrofoam" in item_name:
            best_path, priority = "Mineralization", "Environmental Defense"
            reason = "Waste Trap detected. Recycling is energy-prohibitive. Mineralization ensures zero-microplastic residue."
        elif ar > 90:
            best_path, priority = "Mechanical Recycling", "Resource Preservation"
            reason = "High-value circularity. Resource recovery saves more energy than carbon-return."
        else:
            best_path, priority = "Mineralization", "Environmental Defense"
            reason = "Efficiency gap. Mineralization is the most cost-effective path for this molecular profile."

        return round(br, 1), round(ar, 1), round(bm, 1), round(am, 1), best_path, priority, reason, mw
    except:
        return None

# --- 4. THE INTERFACE ---
st.set_page_config(page_title="VoraCycle Strategic Arbiter", layout="wide")
st.title("🔮 Wraith VoraCycle: Strategic Arbiter")
st.markdown("### *Forensic Audit: Determining the Apex Move for Environment & Enterprise*")

col_search1, col_search2 = st.columns(2)
with col_search1:
    dropdown_query = st.selectbox("🧬 Preselected Problematic Items:", list(product_inventory.keys()))
with col_search2:
    manual_query = st.text_input("✍️ Manual Search (Type any item):")

query = manual_query if manual_query else dropdown_query

if query and query != "Select a problematic item...":
    smiles = product_inventory.get(query, "CCCCCCCC") 
    audit = run_strategic_audit(query, smiles)
    
    if audit:
        br, ar, bm, am, best_path, priority, reason, mw = audit
        
        st.divider()
        st.header(f"🏆 Best Strategic Path: {best_path}")
        st.info(f"**Primary Objective:** {priority} | **Logic:** {reason}")
        
        # DUAL-PATH PERFORMANCE
        st.header("📊 Comparative Forensic Performance")
        p1, p2 = st.columns(2)
        with p1:
            st.subheader("♻️ Path 1: Mechanical Recycling")
            st.metric("After VoraCycle", f"{ar}%", delta=f"Baseline: {br

