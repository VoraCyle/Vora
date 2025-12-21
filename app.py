import streamlit as st
import streamlit_authenticator as st_auth
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
import pandas as pd

# --- 1. Credentials & Auth ---
credentials = {'usernames': {'wraith': {'name': 'Wraith', 'password': 'Vora1630'}}}
authenticator = st_auth.Authenticate(credentials, "vora_cookie", "auth_key", cookie_expiry_days=30)

# --- 2. Login UI ---
authenticator.login(location='main')

# --- 3. Main Application ---
if st.session_state.get("authentication_status"):
    authenticator.logout('Logout', 'sidebar')
    
    st.title("🔮 Wraith VoraCycle")
    st.markdown("### National Wholesaler Strategic Intelligence Dashboard")

    # --- 4. Logic & Grading Engine ---
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
        
        recycle = max(5, min(98, 92 - (rings * 15) - (toxic * 40)))
        fate = max(5, min(98, 20 + (15 if has_bio else 0) - (rings * 10) - (toxic * 30)))
        
        return {"name": name, "mol": mol, "recycle": int(recycle), "fate": int(fate), "toxic": toxic}

    smiles_dict = {
        "Polyethylene (PE) - Bags": "CCCCCCCCCC",
        "PET (Polyester) - Trays": "C1=CC=C(C=C1)C(=O)OCCOC(=O)C2=CC=C(C=C2)C(=O)O",
        "Polypropylene (PP) - Tubs": "CC(C)CCCCCC",
        "Polystyrene (PS) - Foam": "C1=CC=C(C=C1)C(C)C",
        "PVC - Meat Cling Film": "C=CCl",
        "Nylon (PA6) - Vacuum Seals": "CCCCCCN C(=O)CCCCC C(=O)N",
        "PFAS - Grease-proof Paper": "C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F"
    }

    tab1, tab2 = st.tabs(["🔍 Deep Dive Audit", "🌎 Global Market Benchmarking"])

    # --- TAB 1: DEEP DIVE (STABLE & UNTOUCHED) ---
    with tab1:
        st.sidebar.header("Audit Controls")
        category = st.sidebar.selectbox("Application", ["Hot Food", "Cold Storage", "Dry Goods", "Industrial"])
        selected_item = st.selectbox("Select Target Material", list(smiles_dict.keys()) + ["Custom SMILES"])
        
        smiles_input = st.text_input("SMILES Barcode", smiles_dict.get(selected_item, "C1=CC=C(C=C1)C=C"))
        current = analyze_material(smiles_input, selected_item)
        
        if current:
            st.image(Draw.MolToImage(current['mol'], size=(400, 400)), caption=f"Molecular DNA: {selected_item}")
            st.markdown("---")
            st.header("⚖️ The Transformation Analysis")
            red_rec = min(98, current['recycle'] + 25)
            red_fate = min(98, current['fate'] + 45)
            col_b, col_a = st.columns(2)
            with col_b:
                st.markdown("### 🔴 BEFORE (Current)")
                st.metric("Recycle Score", f"{current['recycle']}/100", f"Grade {get_letter_grade(current['recycle'])}", delta_color="inverse")
                st.metric("Landfill Fate", f"{current['fate']}/100", f"Grade {get_letter_grade(current['fate'])}", delta_color="inverse")
            with col_a:
                st.markdown("### 🟢 AFTER (VoraCycle)")
                st.metric("Recycle Score", f"{red_rec}/100", f"Grade {get_letter_grade(red_rec)}")
                st.metric("Landfill Fate", f"{red_fate}/100", f"Grade {get_letter_grade(red_fate)}")
            st.markdown("---")
            st.header("🎯 Procurement Strategy Verdict")
            if category == "Hot Food" or current['recycle'] < 55:
                st.warning("🏁 **STRATEGIC CHOICE: LANDFILL SAFETY (BIO-ASSIMILATION)**")
                st.write("""
                **Why this Verdict?** Items like greasy chicken bags are rejected by recyclers. 
                **Benefit:** VoraCycle redesign ensures **Bio-Mineralization**, meaning the plastic turns into soil rather than microplastics.
                """)
            else:
                st.success("🏁 **STRATEGIC CHOICE: RECYCLE (CIRCULAR RECOVERY)**")
                st.write("""
                **Why this Verdict?** High purity detected. 
                **Benefit:** This material can be sold back into the supply chain, reducing Net Cost of Goods.
                """)

    # --- TAB 2: GLOBAL BENCHMARKING & COMPARISONS ---
    with tab2:
        st.subheader("📊 Global Market Alignment (0-100)")
        if current:
            benchmarks = {
                "Current Costco Item": current['recycle'],
                "Sam's Club Baseline": 58,
                "Walmart Sustainable Goal": 65,
                "EU Grade A standard": 88
            }
            for entity, score in benchmarks.items():
                c_label, c_bar, c_val = st.columns([2, 5, 1])
                c_label.write(f"**{entity}**")
                c_bar.progress(score / 100)
                c_val.write(f"**{score}** ({get_letter_grade(score)})")

        st.markdown("---")
        st.subheader("📋 Multi-Item Outcome Comparison")
        comparison_list = st.multiselect("Benchmark materials side-by-side", list(smiles_dict.keys()), default=list(smiles_dict.keys())[:3])
        
        if comparison_list:
            comp_rows = []
            for item_name in comparison_list:
                res = analyze_material(smiles_dict[item_name], item_name)
                comp_rows.append({
                    "Material": item_name,
                    "Recycle Rating": f"{res['recycle']} ({get_letter_grade(res['recycle'])})",
                    "Landfill Fate": f"{res['fate']} ({get_letter_grade(res['fate'])})",
                    "Safety": "PASS" if res['toxic'] == 0 else "FAIL"
                })
            st.table(pd.DataFrame(comp_rows))

        st.markdown("---")
        st.subheader("📝 Final Strategic Thoughts")
        c1, c2 = st.columns(2)
        with c1:
            st.info("**Why Change?** Avoids 'Plastic Taxes' and eliminates toxins (PFAS) leaching into food fats.")
        with c2:
            st.warning("**Market Gap:** Costco is currently trailing Grade A standards. Redesigning for 'Honest Endgames' is the only way to lead.")

# --- Authentication Logic (Stable) ---
elif st.session_state.get("authentication_status") is False:
    st.error('Login Failed.')
elif st.session_state.get("authentication_status") is None:
    st.warning('Please log in.')
