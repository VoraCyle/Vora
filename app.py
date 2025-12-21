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

    # --- 3. Sidebar: Inventory & Library ---
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

    # Comparison Mode Selection
    mode = st.radio("Select Mode", ["Single Item Diagnostic", "Side-by-Side Comparison"])

    # --- 4. Core Logic Function ---
    def analyze_material(smiles_str, name):
        mol = Chem.MolFromSmiles(smiles_str)
        if not mol:
            return None
        
        # Chemical Descriptors
        mw = Descriptors.MolWt(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        lability = len([a for a in mol.GetAtoms() if a.GetSymbol() in ['O', 'N']])
        toxic = len([a for a in mol.GetAtoms() if a.GetSymbol() in ['Cl', 'F', 'Br']])
        
        # Scoring Engines
        tsi = (rings * 30) + (mw / 5)  # Thermal Stability
        recycle_score = max(5, min(98, 92 - (rings * 12) - (toxic * 35)))
        fate_score = max(5, min(98, 20 + (15 if lability > 0 else 0) - (rings * 10) - (toxic * 25)))
        carbon_savings = 60 if fate_score > 50 else 0 # % vs virgin plastic
        
        return {
            "name": name,
            "mol": mol,
            "mw": mw,
            "tsi": tsi,
            "recycle": recycle_score,
            "fate": fate_score,
            "carbon": carbon_savings,
            "toxic": toxic
        }

    def get_grade(s):
        if s > 80: return "A (Superior)"
        if s > 60: return "B (Standard)"
        if s > 40: return "C (Poor)"
        return "F (Critical)"

    # --- 5. EXECUTION: SINGLE ITEM ---
    if mode == "Single Item Diagnostic":
        selected_name = st.selectbox("Select Item", list(smiles_dict.keys()) + ["Custom SMILES"])
        smiles_input = st.text_input("SMILES String", smiles_dict.get(selected_name, "C1=CC=C(C=C1)C=C"))
        
        data = analyze_material(smiles_input, selected_name)
        if data:
            st.image(Draw.MolToImage(data['mol'], size=(600, 600)))
            
            # Outcome Paths
            st.markdown("---")
            st.header("🏁 The Endgame Strategy")
            p1, p2, p3 = st.columns(3)
            p1.metric("Recycle Path", get_grade(data['recycle']))
            p2.metric("Landfill Path", get_grade(data['fate']))
            p3.metric("CO2 Savings", f"{data['carbon']}%")

            # Safety & Integrity
            st.markdown("---")
            st.subheader("🛡️ Safety & Integrity Audit")
            if category == "Hot Food (Meat/Chicken)" and data['tsi'] < 50:
                st.error("⚠️ **Thermal Risk:** Potential for chemical migration into hot food.")
            else:
                st.success("✅ **Integrity Verified:** Safe for wholesale application temperatures.")
            
            if data['toxic'] > 0:
                st.error(f"🚨 **Toxicity Alert:** Detected {data['toxic']} halogen atoms. High environmental liability.")
            
            # Supplier Audit Button
            if st.button("📄 Generate Supplier Audit Report"):
                report = f"REPORT: {selected_name}\nRECYCLE: {get_grade(data['recycle'])}\nFATE: {get_grade(data['fate'])}\nCARBON SAVED: {data['carbon']}%\nMANDATE: Transition to Bio-Aromatic Furan chains."
                st.text_area("Audit Output:", value=report)

    # --- 6. EXECUTION: SIDE-BY-SIDE ---
    else:
        st.subheader("📊 Comparative Market Benchmarking")
        choices = st.multiselect("Select materials to compare", list(smiles_dict.keys()), default=list(smiles_dict.keys())[:2])
        
        comp_results = []
        for c in choices:
            res = analyze_material(smiles_dict[c], c)
            comp_results.append(res)
        
        if comp_results:
            df = pd.DataFrame(comp_results)
            # Create a clean comparison table
            display_df = df[['name', 'mw', 'recycle', 'fate', 'carbon']].copy()
            display_df.columns = ['Material', 'Weight', 'Recycle Score', 'Landfill Fate', 'CO2 % Saved']
            st.table(display_df)
            
            # Benchmarking against competitors
            st.markdown("---")
            st.subheader("🎯 Competitive Landscape")
            for res in comp_results:
                st.write(f"**{res['name']}** vs Market Average (Walmart/Sam's)")
                st.progress(res['recycle'] / 100)

elif st.session_state.get("authentication_status") is False:
    st.error('Username/password incorrect')
elif st.session_state.get("authentication_status") is None:
    st.warning('Please login to access the VoraCycle Engine.')
