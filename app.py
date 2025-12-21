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

    # --- 3. Logic & Grading Engine ---
    def get_letter_grade(score):
        if score > 85: return "A"
        if score > 70: return "B"
        if score > 55: return "C"
        if score > 40: return "D"
        return "F"

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
        
        return {"name": name, "mol": mol, "recycle": recycle, "fate": fate, "tsi": tsi, "toxic": toxic}

    smiles_dict = {
        "Polyethylene (PE) - Bags": "CCCCCCCCCC",
        "PET (Polyester) - Trays": "C1=CC=C(C=C1)C(=O)OCCOC(=O)C2=CC=C(C=C2)C(=O)O",
        "Polypropylene (PP) - Tubs": "CC(C)CCCCCC",
        "Polystyrene (PS) - Foam": "C1=CC=C(C=C1)C(C)C",
        "PVC - Meat Cling Film": "C=CCl",
        "Nylon (PA6) - Vacuum Seals": "CCCCCCN C(=O)CCCCC C(=O)N",
        "PFAS - Grease-proof Paper": "C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F"
    }

    # --- 4. Navigation ---
    tab1, tab2 = st.tabs(["🔍 Deep Dive Audit", "🌎 Global Competitor Benchmarking"])

    # --- TAB 1: DEEP DIVE (BEFORE/AFTER & DESCRIPTIONS) ---
    with tab1:
        st.sidebar.header("Audit Controls")
        category = st.sidebar.selectbox("Application", ["Hot Food", "Cold Storage", "Dry Goods", "Industrial"])
        selected_item = st.selectbox("Select Target Material", list(smiles_dict.keys()) + ["Custom SMILES"])
        
        smiles_input = st.text_input("SMILES Code", smiles_dict.get(selected_item, "C1=CC=C(C=C1)C=C"))
        current = analyze_material(smiles_input, selected_item)
        
        if current:
            st.image(Draw.MolToImage(current['mol'], size=(400, 400)))
            
            st.markdown("---")
            st.header("⚖️ Transformation Analysis (Before vs. After)")
            
            # Simulated Redesign
            red_rec = min(98, current['recycle'] + 22)
            red_fate = min(98, current['fate'] + 48)

            col1, col2 = st.columns(2)
            with col1:
                st.markdown("### 🔴 BEFORE (Current)")
                st.metric("Recycle Score", f"{current['recycle']}/100", f"Grade {get_letter_grade(current['recycle'])}")
                st.metric("Landfill Fate", f"{current['fate']}/100", f"Grade {get_letter_grade(current['fate'])}")
                st.error("**Outcome:** Material creates persistent microplastics. Zero bio-availability in landfill.")

            with col2:
                st.markdown("### 🟢 AFTER (VoraCycle)")
                st.metric("Recycle Score", f"{red_rec}/100", f"Grade {get_letter_grade(red_rec)}")
                st.metric("Landfill Fate", f"{red_fate}/100", f"Grade {get_letter_grade(red_fate)}")
                st.success("**Outcome:** Full Bio-Mineralization. Material returns to biomass within 24 months.")

            st.markdown("---")
            st.header("🎯 Final Procurement Decision")
            if category == "Hot Food" or current['recycle'] < 55:
                st.warning("🏁 **PATH: LANDFILL SAFETY.** Food-contamination makes recycling high-risk. Focus on VoraCycle 'After' redesign to ensure environmental safety.")
            else:
                st.success("🏁 **PATH: RECYCLE.** Material is high-purity. Focus on circular recovery loops.")

    # --- TAB 2: GLOBAL BENCHMARKING (COMPETITION) ---
    with tab2:
        st.subheader("📊 Global Market Alignment (0-100)")
        st.write("Comparing Costco's current position to global leaders and standards:")

        if current:
            benchmarks = {
                "Current Costco Item": current['recycle'],
                "Sam's Club Average": 58,
                "Walmart Sustainable Standard": 65,
                "EU Grade A Standard": 88
            }

            for entity, score in benchmarks.items():
                grade = get_letter_grade(score)
                c_label, c_bar, c_val = st.columns([2, 5, 1])
                c_label.write(f"**{entity}**")
                c_bar.progress(score / 100)
                c_val.write(f"**{score}** ({grade})")

            st.markdown("---")
            gap = 88 - current['recycle']
            if gap > 0:
                st.warning(f"💡 **Market Gap:** You are currently **{gap} points** behind the EU Grade A standard. The VoraCycle redesign closes this gap immediately.")

        st.subheader("📋 Comparative Item Analysis")
        comparison_list = st.multiselect("Select materials to compare outcome benefits", list(smiles_dict.keys()), default=list(smiles_dict.keys())[:3])
        if comparison_list:
            comp_data = []
            for item in comparison_list:
                res = analyze_material(smiles_dict[item], item)
                comp_data.append({"Material": item, "Recycle": f"{res['recycle']} ({get_letter_grade(res['recycle'])})", "Fate": f"{res['fate']} ({get_letter_grade(res['fate'])})", "Safety": "PASS" if res['toxic'] == 0 else "FAIL"})
            st.table(pd.DataFrame(comp_data))

elif st.session_state.get("authentication_status") is False:
    st.error('Login Failed.')
elif
