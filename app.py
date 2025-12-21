import streamlit as st
import streamlit_authenticator as st_auth
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
import pandas as pd

# --- 1. Credentials & Auth Setup ---
credentials = {'usernames': {'wraith': {'name': 'Wraith', 'password': 'Vora1630'}}}
authenticator = st_auth.Authenticate(credentials, "vora_cookie", "auth_key", cookie_expiry_days=30)

# --- 2. Login UI ---
authenticator.login(location='main')

if st.session_state.get("authentication_status"):
    authenticator.logout('Logout', 'sidebar')
    
    st.title("🔮 Wraith VoraCycle")
    st.markdown("### National Wholesaler Diagnostic & Procurement Engine")

    # --- 3. Sidebar Configuration ---
    st.sidebar.header("🏢 Warehouse Inventory")
    category = st.sidebar.selectbox(
        "Application Type",
        ["Hot Food (Meat/Chicken)", "Cold Storage (Produce/Dairy)", "Dry Goods (Pantry)", "Logistics (Pallets/Stretch)", "Industrial (Cleaning)"]
    )

    smiles_dict = {
        "Polyethylene (PE) - Bags": "CCCCCCCCCC",
        "PET (Polyester) - Trays": "C1=CC=C(C=C1)C(=O)OCCOC(=O)C2=CC=C(C=C2)C(=O)O",
        "Polypropylene (PP) - Tubs": "CC(C)CCCCCC",
        "Polystyrene (PS) - Foam": "C1=CC=C(C=C1)C(C)C",
        "PVC - Cling Film": "C=CCl",
        "Nylon (PA6) - Meat Seals": "CCCCCCN C(=O)CCCCC C(=O)N",
        "PFAS - Grease-proof Paper": "C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F"
    }

    mode = st.radio("Select Mode", ["Single Item Diagnostic", "Side-by-Side Comparison"])

    # --- 4. Logic Functions ---
    def get_grade(s):
        if s > 80: return "A (Superior)"
        if s > 60: return "B (Standard)"
        if s > 40: return "C (Poor)"
        return "F (Critical)"

    def analyze_material(smiles_str, name):
        mol = Chem.MolFromSmiles(smiles_str)
        if not mol: return None
        
        mw = Descriptors.MolWt(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        lability = len([a for a in mol.GetAtoms() if a.GetSymbol() in ['O', 'N']])
        toxic = len([a for a in mol.GetAtoms() if a.GetSymbol() in ['Cl', 'F', 'Br']])
        
        tsi = (rings * 30) + (mw / 5)  # Thermal Stability
        recycle_score = max(5, min(98, 92 - (rings * 12) - (toxic * 35)))
        fate_score = max(5, min(98, 20 + (15 if lability > 0 else 0) - (rings * 10) - (toxic * 25)))
        carbon_savings = 60 if fate_score > 50 else 0 
        
        return {
            "name": name, "mol": mol, "mw": mw, "tsi": tsi, 
            "recycle": recycle_score, "fate": fate_score, 
            "carbon": carbon_savings, "toxic": toxic
        }

    def get_rationale(data, category):
        reasons = []
        if data['toxic'] > 0:
            reasons.append("❌ **BAD IDEA:** Current material contains halogens. These can migrate into fats (leaching) and release toxic gases if incinerated.")
        if category == "Hot Food (Meat/Chicken)" and data['tsi'] < 50:
            reasons.append("❌ **BAD IDEA:** Low thermal resistance. The packaging may soften or 'bond' to the food at 180°F, contaminating the flavor and safety.")
        
        if data['fate'] < 50:
            reasons.append("💡 **THE CHANGE:** Transition to a 'Bio-Aromatic' backbone (like PEF).")
            reasons.append("🛡️ **THE BENEFIT:** This matches the strength of plastic for heavy items but allows soil microbes to 'unzip' the material in a landfill.")
        
        if data['recycle'] > 75:
            reasons.append("✅ **GOOD IDEA:** High-purity material. Reduces Plastic Tax liability.")
        return reasons

    # --- 5. Main Display Logic ---
    if mode == "Single Item Diagnostic":
        selected_name
