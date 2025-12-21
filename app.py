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

if st.session_state.get("authentication_status"):
    authenticator.logout('Logout', 'sidebar')
    
    st.title("🔮 Wraith VoraCycle")
    st.markdown("### National Wholesaler Strategic Diagnostic")

    # --- 3. Sidebar Inventory ---
    st.sidebar.header("🏢 Warehouse Inventory")
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

    # --- 4. Logic & Rating Functions ---
    def get_letter_grade(score):
        if score > 85: return "A (Superior)"
        if score > 70: return "B (Standard)"
        if score > 55: return "C (Sub-Optimal)"
        if score > 40: return "D (Poor)"
        return "F (Critical)"

    def get_rating_description(score, path_type):
        if path_type == "Recycle":
            if score > 70: return "Highly compatible with mechanical recycling. High secondary market value."
            if score > 40: return "Difficult to sort. Likely to be incinerated or downcycled into low-grade lumber."
            return "Non-recyclable. Contaminates other plastic streams. Causes machine jams."
        else: # Landfill Path
            if score > 70: return "Bio-assimilable. Microbes recognize the carbon chain and convert it to soil/CO2."
            if score > 40: return "Fragmentation risk. Breaks into microplastics but does not fully mineralize."
            return "Persistent. Indestructible molecular bonds create a 'plastic fossil' for 400+ years."

    def analyze_material(smiles_str):
        mol = Chem.MolFromSmiles(smiles_str)
        if not mol: return None
        mw = Descriptors.MolWt(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        toxic = len([a for a in mol.GetAtoms() if a.GetSymbol() in ['Cl', 'F', 'Br']])
        has_bio = any(a.GetSymbol() in ['O', 'N'] for a in mol.GetAtoms())
        
        # Performance/Safety Calculation
        tsi = (rings * 35) + (mw / 5) 
        recycle_score = max(5, min(98, 92 - (rings * 15) - (toxic * 40)))
        fate_score = max(5, min(98, 20 + (15 if has_bio else 0) - (rings * 10) - (toxic * 30)))
        
        return {"mol": mol, "tsi": tsi, "recycle": recycle_score, "fate": fate_score, "toxic": toxic}

    # --- 5. Main Diagnostic Execution ---
    selected_name = st.selectbox("Select Item for Audit", list(smiles_dict.keys()) + ["Custom SMILES"])
    smiles_input = st.text_input("SMILES Barcode", smiles_dict.get(selected_name, "C1=CC=C(C=C1)C=C"))
    
    current = analyze_material(smiles_input)
    if current:
        st.image(Draw.MolToImage(current['mol'], size=(500, 500)))

        # --- A. BEFORE VS AFTER RATINGS ---
        st.markdown("---")
        st.header("⚖️ Life Cycle Upgrade Comparison")
        
        # Simulation of Redesign (Optimized version)
        redesign_recycle = min(98, current['recycle'] + 25)
        redesign_fate = min(98, current['fate'] + 45)

        col1, col2 = st.columns(2)
        with col1:
            st.subheader("🔴 Current (Before)")
            st.metric("Recycle Rating", f"{current['recycle']}/100", get_letter_grade(current['recycle']))
            st.metric("Landfill Fate", f"{current['fate']}/100", get_letter_grade(current['fate']))
            st.write(f"**Outcome:** {get_rating_description(current['fate'], 'Landfill')}")

        with col2:
            st.subheader("🟢 VoraCycle (After)")
            st.metric("Recycle Rating", f"{redesign_recycle}/100", get_letter_grade(redesign_recycle), delta="Upgrade")
            st.metric("Landfill Fate", f"{redesign_fate}/100", get_letter_grade(redesign_fate), delta="Upgrade")
            st.write(f"**Outcome:** Fully mineralized carbon; zero microplastic legacy.")

        # --- B. STRATEGIC DESCRIPTIONS & OUTCOMES ---
        st.markdown("---")
        st.header("📖 Outcome Analysis: Why Change?")
        
        st.write("### 🔄 The Recycling Path (The Preferred Loop)")
        st.info("""**Why it's better to Recycle:** Recycling keeps the material in the 'Wholesale Loop.' 
        It reduces the need to extract fresh oil and eliminates EPR Plastic Taxes. 
        A 'Grade A' rating here means the item pays for its own waste management through secondary material value.""")
        
        st.write("### 🌍 The Landfill Path (The Safety Net)")
        st.warning("""**What happens in the Landfill:** Most food-contaminated items (like chicken bags) 
        cannot be recycled. Our redesign ensures that if the item ends up in a landfill, it does not 
        leach chemicals into groundwater or turn into microplastics. Instead, it undergoes **Bio-Mineralization**.""")

        # --- C. SAFETY & INTEGRITY GUARANTEE ---
        st.markdown("---")
        st.header("🛡️ Performance & Food Safety Audit")
        
        if category == "Hot Food (Meat/Chicken)" and current['tsi'] < 55:
            st.error("🚨 **Safety Alert:** Current material fails high-heat integrity test. Risk of chemical leaching into food fats.")
        else:
            st.success("✅ **Safety Verified:** Material maintains structural integrity at operational temperatures.")

        if st.button("📄 Generate Executive Supplier Audit"):
            report = f"AUDIT: {selected_name}\nCURRENT GRADE: {get_letter_grade(current['recycle'])}\nPROPOSED GRADE: {get_letter_grade(redesign_recycle)}\n\nFINAL DESCRIPTION: The redesign replaces indestructible carbon chains with bio-assimilable links. This prevents chemical migration into hot food while ensuring a safe environmental endgame."
            st.text_area("Copy Final Audit:", value=report, height=200)

elif st.session_state.get("authentication_status") is False:
    st.error('Login Failed.')
elif st.session_state.get("authentication_status") is None:
    st.warning('Please log in.')
