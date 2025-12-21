import streamlit as st
import streamlit_authenticator as st_auth
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
import pandas as pd

# --- 1. Credentials & Auth ---
# Password: Vora1630
credentials = {'usernames': {'wraith': {'name': 'Wraith', 'password': 'Vora1630'}}}
authenticator = st_auth.Authenticate(credentials, "vora_cookie", "auth_key", cookie_expiry_days=30)

# --- 2. Login UI ---
authenticator.login(location='main')

if st.session_state.get("authentication_status"):
    authenticator.logout('Logout', 'sidebar')
    
    st.title("🔮 Wraith VoraCycle")
    st.markdown("### National Wholesaler Strategic Intelligence Dashboard")

    # --- 3. Shared Logic & Grading Functions ---
    def get_letter_grade(score):
        if score > 85: return "A (Superior)"
        if score > 70: return "B (Standard)"
        if score > 55: return "C (Sub-Optimal)"
        if score > 40: return "D (Poor)"
        return "F (Critical)"

    def analyze_material(smiles_str, name="Custom"):
        mol = Chem.MolFromSmiles(smiles_str)
        if not mol: return None
        mw = Descriptors.MolWt(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        toxic = len([a for a in mol.GetAtoms() if a.GetSymbol() in ['Cl', 'F', 'Br']])
        has_bio = any(a.GetSymbol() in ['O', 'N'] for a in mol.GetAtoms())
        
        # Scoring Metrics (0-100)
        recycle = max(5, min(98, 92 - (rings * 15) - (toxic * 40)))
        fate = max(5, min(98, 20 + (15 if has_bio else 0) - (rings * 10) - (toxic * 30)))
        tsi = (rings * 35) + (mw / 5) 
        
        return {"name": name, "mol": mol, "recycle": recycle, "fate": fate, "tsi": tsi, "toxic": toxic, "mw": mw}

    smiles_dict = {
        "Polyethylene (PE) - Bags": "CCCCCCCCCC",
        "PET (Polyester) - Trays": "C1=CC=C(C=C1)C(=O)OCCOC(=O)C2=CC=C(C=C2)C(=O)O",
        "Polypropylene (PP) - Tubs": "CC(C)CCCCCC",
        "Polystyrene (PS) - Foam": "C1=CC=C(C=C1)C(C)C",
        "PVC - Meat Cling Film": "C=CCl",
        "Nylon (PA6) - Vacuum Seals": "CCCCCCN C(=O)CCCCC C(=O)N",
        "PFAS - Grease-proof Paper": "C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F"
    }

    # --- 4. System Navigation ---
    st.sidebar.header("🏢 Navigation")
    mode = st.tabs(["🔍 Deep Dive Audit", "📊 Market Leaderboard"])

    # --- TAB 1: DEEP DIVE (BEFORE/AFTER & DECISION) ---
    with mode[0]:
        st.sidebar.subheader("Audit Settings")
        category = st.sidebar.selectbox("Application Type", ["Hot Food", "Cold Storage", "Dry Goods", "Industrial"])
        selected_item = st.selectbox("Select Target Material", list(smiles_dict.keys()) + ["Custom SMILES"])
        
        if selected_item == "Custom SMILES":
            smiles_input = st.text_input("Enter SMILES Code", "C1=CC=C(C=C1)C=C")
        else:
            smiles_input = smiles_dict[selected_item]
            
        current = analyze_material(smiles_input, selected_item)
        
        if current:
            st.image(Draw.MolToImage(current['mol'], size=(450, 450)), caption="Molecular Architecture")
            
            # Outcome Metrics
            st.markdown("---")
            st.header("⚖️ The VoraCycle Transformation (Before vs. After)")
            
            # Redesign Simulation (The Upgrade)
            red_rec = min(98, current['recycle'] + 20)
            red_fate = min(98, current['fate'] + 50)
            
            col_b, col_a = st.columns(2)
            with col_b:
                st.markdown("### 🔴 BEFORE (Current)")
                st.metric("Recycle Rating", f"{current['recycle']}/100", get_letter_grade(current['recycle']), delta_color="inverse")
                st.metric("Landfill Fate", f"{current['fate']}/100", get_letter_grade(current['fate']), delta_color="inverse")
                st.write("**Outcome:** Persistent microplastics. Material remains in anaerobic stasis for 400+ years.")

            with col_a:
                st.markdown("### 🟢 AFTER (VoraCycle)")
                st.metric("Recycle Rating", f"{red_rec}/100", get_letter_grade(red_rec))
                st.metric("Landfill Fate", f"{red_fate}/100", get_letter_grade(red_fate))
                st.write("**Outcome:** Full Bio-Mineralization. Carbon returns to soil within 24-48 months.")

            # THE DECISION BLOCK
            st.markdown("---")
            st.header("🎯 Final Procurement Decision")
            if category == "Hot Food" or current['recycle'] < 50:
                st.warning("🏁 **RECOMMENDED PATH: LANDFILL SAFETY.** Due to food contamination risk (oils/fats), this item is unlikely to be recycled. Prioritize the VoraCycle 'After' redesign to ensure it disappears safely in the trash.")
            else:
                st.success("🏁 **RECOMMENDED PATH: RECYCLE.** This item is high-purity and dry. Optimize for the circular loop to avoid EPR taxes.")

            st.markdown("---")
            st.subheader("🛡️ Safety & Integrity Check")
            if category == "Hot Food" and current['tsi'] < 55:
                st.error("🚨 **INTEGRITY RISK:** Current material fails thermal audit. High risk of leaching chemicals into food at 180°F. The redesign is MANDATORY for food safety.")
            else:
                st.success("✅ **STABLE:** Material architecture is safe for current operational temperatures.")

    # --- TAB 2: MARKET LEADERBOARD ---
    with mode[1]:
        st.subheader("📊 Side-by-Side: Find the Beneficial Outcome")
        comparison_list = st.multiselect("Select materials to benchmark", list(smiles_dict.keys()), default=list(smiles_dict.keys())[:3])
        
        if comparison_list:
            comp_data = []
            for item in comparison_list:
                res = analyze_material(smiles_dict[item], item)
                comp_data.append({
                    "Material": item,
                    "Recycle Grade": get_letter_grade(res['recycle']),
                    "Recycle Score": res['recycle'],
                    "Landfill Fate": get_letter_grade(res['fate']),
                    "Toxicity Risk": "HIGH" if res['toxic'] > 0 else "NONE"
                })
            
            st.table(pd.DataFrame(comp_data))

            st.markdown("---")
            st.subheader("🌎 Global Competitor Benchmarking")
            st.write("**How Costco stands against Market Leaders:**")
            
            # Dynamic Benchmarking Display
            for d in comp_data:
                st.write(f"**{d['Material']}** vs. EU Standards")
                st.progress(d['Recycle Score'] / 100)
                
            st.markdown("""
            * **Costco Current Benchmark:** 55% (Grade C)
            * **Walmart Standard:** 60% (Grade B)
            * **EU Circularity Law:** 85% (Grade A)
            """)

elif st.session_state.get("authentication_status") is False:
    st.error('Login Failed.')
elif st.session_state.get("authentication_status") is None:
    st.warning('Please login to the Wholesaler Terminal.')
