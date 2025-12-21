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

    # --- TAB 1: DEEP DIVE & STRATEGIC RATIONALE ---
    with tab1:
        st.sidebar.header("Audit Controls")
        category = st.sidebar.selectbox("Application", ["Hot Food", "Cold Storage", "Dry Goods", "Industrial"])
        selected_item = st.selectbox("Select Target Material", list(smiles_dict.keys()) + ["Custom SMILES"])
        
        smiles_input = st.text_input("SMILES Barcode", smiles_dict.get(selected_item, "C1=CC=C(C=C1)C=C"))
        current = analyze_material(smiles_input, selected_item)
        
        if current:
            st.image(Draw.MolToImage(current['mol'], size=(400, 400)), caption=f"Molecular DNA: {selected_item}")
            
            st.markdown("---")
            st.header("⚖️ The Transformation Analysis (Before vs. After)")
            
            # Simulated Redesign Metrics
            red_rec = min(98, current['recycle'] + 25)
            red_fate = min(98, current['fate'] + 45)

            col_b, col_a = st.columns(2)
            with col_b:
                st.markdown("### 🔴 BEFORE (Current)")
                st.metric("Recycle Rating", f"{current['recycle']}/100", f"Grade {get_letter_grade(current['recycle'])}", delta_color="inverse")
                st.metric("Landfill Fate", f"{current['fate']}/100", f"Grade {get_letter_grade(current['fate'])}", delta_color="inverse")
                st.error("**Persistence Liability:** This material creates microplastics that fragment but do not disappear for 400+ years.")

            with col_a:
                st.markdown("### 🟢 AFTER (VoraCycle)")
                st.metric("Recycle Rating", f"{red_rec}/100", f"Grade {get_letter_grade(red_rec)}")
                st.metric("Landfill Fate", f"{red_fate}/100", f"Grade {get_letter_grade(red_fate)}")
                st.success("**Bio-Assimilation:** Carbon returns to the soil safely via mineralization within 24-48 months.")

            st.markdown("---")
            st.header("🎯 Procurement Strategy Verdict")
            if category == "Hot Food" or current['recycle'] < 55:
                st.warning("🏆 **DECISION: LANDFILL SAFETY.** Food-grease contamination makes recycling unviable. Transitioning to a Bio-Aromatic polymer ensures a safe environmental return.")
            else:
                st.success("🏆 **DECISION: RECYCLE.** High purity detected. Focus on circular recovery and 'Closed-Loop' logistics.")

    # --- TAB 2: GLOBAL BENCHMARKING & FINAL SUGGESTIONS ---
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
            st.subheader("📝 Final Strategic Suggestions")
            
            st.write("#### 1. Why Change is Beneficial in All Aspects")
            st.info("""
            * **Financial:** Avoids upcoming 'Plastic Taxes' (approx. $420/tonne for Grade F materials).
            * **Safety:** Eliminates halogenated toxins (PFAS/PVC) that leach into rotisserie chicken and meat fats.
            * **Brand:** Moves Costco from 'Compliance' to 'Leadership,' justifying membership fees through superior Kirkland standards.
            """)

            st.write("#### 2. Competitor Benchmarking Suggestions")
            st.warning(f"""
            * **The Gap:** You are currently **{88 - current['recycle']} points** behind the EU Grade A standard.
            * **Action:** Use the VoraCycle 'After' redesign to leapfrog domestic competitors (Walmart/Sam's) and hit the 85+ mark. This standardizes your supply chain globally and prevents local regulatory shocks.
            """)

# --- 5. Authentication Error Handling (Fixed) ---
elif st.session_state.get("authentication_status") is False:
    st.error('Login Failed. Please check your credentials.')

elif st.session_state.get("authentication_status") is None:
    st.warning('Please log in to the VoraCycle Dashboard.')
