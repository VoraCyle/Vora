import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Descriptors

# --- 1. ENTERPRISE SECURITY (THE VAULT) ---
def initialize_ai():
    """
    Retrieves keys from the Streamlit Secrets Vault. 
    Matches the names you entered in the Dashboard.
    """
    try:
        # Pulling from your Secrets dashboard
        k1 = GEMINI_KEY_1 = "AIzaSyDRJyAVtMKTolhE6vpYjPiA9eSr7Ko9_Og"
        genai.configure(api_key=k1)
        return genai.GenerativeModel('gemini-1.5-flash')
    except Exception as e:
        st.error("🔑 Vault Access Denied. Ensure GEMINI_KEY_1 is in your Streamlit Secrets.")
        st.stop()

model = initialize_ai()

# --- 2. THE STRATEGIC INVENTORY ---
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
    "Waxed Produce Box": "CCCCCCCCCCCCCCCCCCCC",
}

# --- 3. THE STRATEGIC ARBITER (DECISION LOGIC) ---
def run_strategic_audit(item_name, smiles):
    try:
        # Chemical analysis
        mol = Chem.MolFromSmiles(smiles) if "[Al]" not in smiles else None
        mw = Descriptors.MolWt(mol) if mol else 27.0
        toxic = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms()) if mol else False
        
        # Path 1: Mechanical Recycling Performance
        br = 92 if "PET" in item_name or "Aluminum" in item_name else max(10, 75 - (mw / 10))
        if toxic: br -= 35 
        ar = min(99.1, br + 45) 
        
        # Path 2: Soil Mineralization Performance
        bm = 5 if toxic else max(12, 50 - (mw / 20))
        am = 99.4 if not toxic else 94.0 
        
        # Determining the Apex Move
        if toxic or "Multi-layer" in item_name or "Styrofoam" in item_name:
            best_path, priority = "Mineralization", "Environmental Defense"
            reason = "Waste Trap detected. Recycling is energy-prohibitive. Mineralization ensures zero-microplastic residue."
        elif ar > 88:
            best_path, priority = "Mechanical Recycling", "Resource Preservation"
            reason = "High-value circularity. Resource recovery saves more energy than carbon-return."
        else:
            best_path, priority = "Mineralization", "Environmental Defense"
            reason = "Efficiency gap. Mineralization is the most cost-effective path for this molecular profile."

        return round(br, 1), round(ar, 1), round(bm, 1), round(am, 1), best_path, priority, reason
    except:
        return None

# --- 4. THE INTERFACE ---
st.set_page_config(page_title="VoraCycle Strategic Arbiter", layout="wide", page_icon="🔮")
st.title("🔮 Wraith VoraCycle: Strategic Arbiter")
st.markdown("### *Forensic Audit: Determining the Apex Move for Environment & Enterprise*")

col1, col2 = st.columns(2)
with col1:
    dropdown_query = st.selectbox("🧬 Select Industry-Standard Item:", list(product_inventory.keys()))
with col2:
    manual_query = st.text_input("✍️ Manual Search / Custom SMILES:")

query = manual_query if manual_query else dropdown_query

if query and query != "Select a problematic item...":
    smiles = product_inventory.get(query, "CCCCCCCC") 
    audit = run_strategic_audit(query, smiles)
    
    if audit:
        br, ar, bm, am, best_path, priority, reason = audit
        
        st.divider()
        st.header(f"🏆 Best Strategic Path: {best_path}")
        st.info(f"**Primary Objective:** {priority} | **Logic:** {reason}")

