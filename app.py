import streamlit as st
import streamlit_authenticator as st_auth
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
import pandas as pd

# --- 1. Credentials & Auth ---
# Password: Vora1630 (Wraith's secure access)
credentials = {'usernames': {'wraith': {'name': 'Wraith', 'password': 'Vora1630'}}}
authenticator = st_auth.Authenticate(credentials, "vora_cookie", "auth_key", cookie_expiry_days=30)

# --- 2. Login UI ---
authenticator.login(location='main')

if st.session_state.get("authentication_status"):
    authenticator.logout('Logout', 'sidebar')
    
    st.title("🔮 Wraith VoraCycle")
    st.markdown("### The 'Start-Line' Diagnostic: Wholesaler Endgame Engine")

    # --- 3. Global Material Library & Benchmarks ---
    st.sidebar.header("🏢 Inventory Controls")
    category = st.sidebar.selectbox(
        "Application Type",
        ["Hot Food (Meat/Chicken)", "Cold Storage (Produce/Dairy)", "Dry Goods (Pantry)", "Industrial/Cleaning"]
    )

    smiles_dict = {
        "Polyethylene (PE) - Bags": "CCCCCCCCCC",
        "PET (Polyester) - Trays": "C1=CC=C(C=C1)C(=O)OCCOC(=O)C2=CC=C(C=C2)C(=O)O",
        "Polypropylene (PP) - Tubs": "CC(C)CCCCCC",
        "Polystyrene (PS) - Foam": "C1=CC=C(C=C1)C(C)C",
        "PVC - Meat Cling Film": "C=CCl",
        "Nylon (PA6) - Vacuum Seals": "CCCCCCN C(=O)CCCCC C(=O)N",
        "PFAS - Grease-proof Paper": "C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F"
    }

    mode = st.radio("System Mode", ["Deep Diagnostic", "Side-by-Side Market Benchmarking"])

    # --- 4. The Engineering Functions ---
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
        toxic = len([a for a in mol.GetAtoms() if a.GetSymbol() in ['Cl', 'F', 'Br']])
        has_bio_handle = any(a.GetSymbol() in ['O', 'N'] for a in mol.GetAtoms())
        
        tsi = (rings * 35) + (mw / 5)  # Thermal Stability Index
        recycle_score = max(5, min(98, 92 - (rings * 15) - (toxic * 40)))
        fate_score = max(5, min(98, 20 + (15 if has_bio_handle else 0) - (rings * 10) - (toxic * 30)))
        
        # 2025 Financial Logic (EPR Fees based on $420/tonne for non-recyclable)
        tax_cost = 0 if recycle_score > 60 else 420 
        
        return {
            "name": name, "mol": mol, "tsi": tsi, 
            "recycle": recycle_score, "fate": fate_score, 
            "toxic": toxic, "tax": tax_cost, "has_bio": has_bio_handle
        }

    def get_rationale(data, category):
        # The 'Logic Layer' for Costco Executives
        notes = []
        if data['toxic'] > 0:
            notes.append("🚨 **CRITICAL RISK:** Contains halogens (PFAS/Chlorine). These leach into food fats and create permanent toxic environmental liability.")
        if category == "Hot Food (Meat/Chicken)" and data['tsi'] < 55:
            notes.append("🚨 **INTEGRITY FAILURE:** Low heat resistance. The material will likely lose structural strength or migrate chemicals into food at 180°F.")
        
        if data['fate'] < 50:
            notes.append("💡 **THE CHANGE:** Transition to a 'Labile' backbone (Ester or Amide links).")
            notes.append("🛡️ **THE BENEFIT:** Allows microbes to digest the bag in a landfill while keeping the product 100% fresh on the shelf.")
        
        if data['recycle'] > 75:
            notes.append("✅ **ADVANTAGE:** High-purity polymer. Eligible for 'Eco-Modulation' tax breaks ($0/tonne EPR fees).")
        return notes

    # --- 5. EXECUTION: SINGLE DIAGNOSTIC ---
    if mode == "Deep Diagnostic":
        selected_name = st.selectbox("Select Material Type", list(smiles_dict.keys()) + ["Custom SMILES"])
        smiles_input = st.text_input("SMILES Barcode", smiles_dict.get(selected_name, "C1=CC=C(C=C1)C=C"))
        
        data = analyze_material(smiles_input, selected_name)
        if data:
            col1, col2 = st.columns([1, 2])
            with col1:
                st.image(Draw.MolToImage(data['mol'], size=(400, 400)))
            with col2:
                st.subheader("📝 Strategic Rationale")
                for n in get_rationale(data, category):
                    st.write(n)

            st.markdown("---")
            st.header("🏁 The Endgame Outcome")
            p1, p2, p3 = st.columns(3)
            p1.metric("Recycle Path", get_grade(data['recycle']), help="How likely a facility is to accept this.")
            p2.metric("Landfill Path", get_grade(data['fate']), help="Reaction in the trash.")
            p3.metric("EPR Tax Liability", f"${data['tax']}/tonne", delta="Avoidable Cost")

            # Reaction Description
            if data['fate'] < 40:
                st.error("**Landfill Reaction:** Anaerobic Stasis. This item will not break down for 400+ years.")
            else:
                st.success("**Landfill Reaction:** Bio-Mineralization. Carbon returns to biomass safely.")

            if st.button("📄 Generate Supplier Audit"):
                report = f"VORACYCLE AUDIT: {selected_name}\nRECYCLE: {get_grade(data['recycle'])}\nLANDFILL: {get_grade(data['fate'])}\nSAFETY: {'FAIL' if data['toxic']>0 else 'PASS'}\nACTION: Redesign to remove Halogenated chains."
                st.text_area("Official Report:", value=report, height=200)

    # --- 6. EXECUTION: SIDE-BY-SIDE ---
    else:
        st.subheader("📊 Competitive Dashboard (Costco vs. Market)")
        choices = st.multiselect("Select 3 Items to Compare", list(smiles_dict.keys()), default=list(smiles_dict.keys())[:3])
        
        res_list = [analyze_material(smiles_dict[c], c) for c in choices]
        if res_list:
            df = pd.DataFrame(res_list)
            st.table(df[['name', 'recycle', 'fate', 'tax']])
            
            # Global Progress Bars
            for r in res_list:
                st.write(f"**{r['name']}** vs. Global Best-in-Class (EU Standard)")
                st.progress(r['recycle'] / 100)

elif st.session_state.get("authentication_status") is False:
    st.error('Access Denied.')
elif st.session_state.get("authentication_status") is None:
    st.warning('Please log in to the VoraCycle Wholesaler Terminal.')
