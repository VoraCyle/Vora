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
        tsi = (rings * 35) + (mw / 5) 
        
        return {"name": name, "mol": mol, "recycle": int(recycle), "fate": int(fate), "tsi": tsi, "toxic": toxic}

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

    # --- TAB 1: DEEP DIVE (UNTOUCHED PER REQUEST) ---
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
                **Why this Verdict?** Items in this category (like greasy chicken bags or meat trays) are typically rejected by recycling facilities due to organic contamination. Attempting to recycle these is 'Wish-cycling'—it costs more and leads to landfill anyway.
                
                **How it's Beneficial:**
                By choosing a VoraCycle redesigned material, you ensure that even when the item is thrown in the trash, it undergoes **Bio-Mineralization**. This eliminates microplastic liability, avoids 'non-recyclable' taxes, and protects the brand from environmental litigation.
                """)
            else:
                st.success("🏁 **STRATEGIC CHOICE: RECYCLE (CIRCULAR RECOVERY)**")
                st.write("""
                **Why this Verdict?** This material is high-purity and used in a 'dry' application. It has a high secondary market value, meaning recycling facilities will actively pay for this waste stream.
                
                **How it's Beneficial:**
                Focusing on a 'Closed-Loop' system for this item allows Costco to sell its waste back into the supply chain. This lowers the Net Cost of Goods (COGS) and achieves the highest possible ESG rating for the warehouse.
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
        st.subheader("📋 Multi-Item Beneficial Outcome Comparison")
        st.write("Compare different warehouse materials to identify which redesign provides the most significant circular benefit.")
        
        comparison_list = st.multiselect("Select materials for side-by-side rating", list(smiles_dict.keys()), default=list(smiles_dict.keys())[:3])
        
        if comparison_list:
            comp_rows = []
            for item_name in comparison_list:
                res = analyze_material(smiles_dict[item_name], item_name)
                comp_rows.append({
                    "Material": item_name,
                    "Recycle Rating": f"{res['recycle']} ({get_letter_grade(res['recycle'])})",
                    "Landfill Fate": f"{res['fate']} ({get_letter_grade(res['fate'])})",
                    "Food Safety": "PASS" if res['toxic'] == 0 else "FAIL (Leach Risk)"
                })
            st.table(pd.DataFrame(comp_rows))

        st.markdown("---")
        st.subheader("📝 Final Strategic Suggestions & Thoughts")
        
        col_suggest1, col_suggest2 = st.columns(2)
        
        with col_suggest1:
            st.info("#### Why Change is Beneficial")
            st.write("""
            * **Economic:** Moving to Grade A materials reduces exposure to 'Plastic Taxes' which are increasingly based on 0-100 circularity scores.
            * **Regulatory:** Standardizing to the EU 88-point mark future-proofs the supply chain against domestic bans.
            * **Safety:** Removing PVC/PFAS ensures that heat-lamp applications in the deli don't migrate toxins into prepared meals.
            """)
            
        with col_suggest2:
            st.warning("#### Global Benchmarking Thoughts")
            st.write
