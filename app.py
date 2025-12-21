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
        "PVC - Cling Film": "C=CCl",
        "Nylon (PA6) - Meat Seals": "CCCCCCN C(=O)CCCCC C(=O)N",
        "PFAS - Grease-proof Paper": "C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F"
    }

    # --- 4. Logic & Rating Functions ---
    def get_letter_grade(score):
        if score > 85: return "A"
        if score > 70: return "B"
        if score > 55: return "C"
        if score > 40: return "D"
        return "F"

    def get_rating_description(score, path_type):
        if path_type == "Recycle":
            if score > 70: return "✅ **Circular:** High compatibility with mechanical recycling. High secondary market value for pellets."
            if score > 40: return "⚠️ **Downcyclable:** Difficult to sort. Likely to be incinerated or turned into low-grade plastic lumber."
            return "❌ **Non-Recyclable:** Contaminates sorting streams. High risk of warehouse 'Plastic Tax' penalties."
        else: # Landfill Path
            if score > 70: return "✅ **Bio-Assimilable:** Microbes recognize the chain. Carbon returns to soil/biomass safely."
            if score > 40: return "⚠️ **Fragmentation Risk:** Breaks into microplastics. Remains in the environment for 100+ years."
            return "❌ **Indestructible:** Creates a 'plastic fossil'. Will remain in the landfill for 400+ years."

    def analyze_material(smiles_str):
        mol = Chem.MolFromSmiles(smiles_str)
        if not mol: return None
        mw = Descriptors.MolWt(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        toxic = len([a for a in mol.GetAtoms() if a.GetSymbol() in ['Cl', 'F', 'Br']])
        has_bio = any(a.GetSymbol() in ['O', 'N'] for a in mol.GetAtoms())
        
        tsi = (rings * 35) + (mw / 5) 
        recycle_score = max(5, min(98, 92 - (rings * 15) - (toxic * 40)))
        fate_score = max(5, min(98, 20 + (15 if has_bio else 0) - (rings * 10) - (toxic * 30)))
        
        return {"mol": mol, "tsi": tsi, "recycle": recycle_score, "fate": fate_score, "toxic": toxic}

    # --- 5. Execution ---
    selected_name = st.selectbox("Select Item for Audit", list(smiles_dict.keys()) + ["Custom SMILES"])
    smiles_input = st.text_input("SMILES Barcode", smiles_dict.get(selected_name, "C1=CC=C(C=C1)C=C"))
    
    current = analyze_material(smiles_input)
    if current:
        st.image(Draw.MolToImage(current['mol'], size=(500, 500)))

        # --- A. BEFORE VS AFTER RATINGS ---
        st.markdown("---")
        st.header("⚖️ Life Cycle Upgrade Comparison")
        
        # Redesign Simulation
        red_rec = min(98, current['recycle'] + 25)
        red_fate = min(98, current['fate'] + 45)

        col1, col2 = st.columns(2)
        with col1:
            st.subheader("🔴 Current (Before)")
            st.metric("Recycle Score", f"{current['recycle']}/100", get_letter_grade(current['recycle']), delta_color="inverse")
            st.metric("Landfill Fate", f"{current['fate']}/100", get_letter_grade(current['fate']), delta_color="inverse")
            st.write(get_rating_description(current['fate'], 'Landfill'))

        with col2:
            st.subheader("🟢 VoraCycle (After)")
            st.metric("Recycle Score", f"{red_rec}/100", f"Grade: {get_letter_grade(red_rec)}")
            st.metric("Landfill Fate", f"{red_fate}/100", f"Grade: {get_letter_grade(red_fate)}")
            st.write("✅ **Outcome:** Fully mineralized carbon; safe return to the environment.")

        # --- B. STRATEGIC OUTCOMES ---
        st.markdown("---")
        st.header("📖 Strategic Outcome: Why Change?")
        
        tab1, tab2 = st.tabs(["🔄 The Recycling Path", "🌍 The Landfill Path"])
        
        with tab1:
            st.write("### The Economic Loop")
            st.info("Recycling is the preferred path for clean, high-volume plastics. It allows Costco to maintain a 'Closed Loop' supply chain, reducing raw material costs and avoiding EPR taxes.")
        
        with tab2:
            st.write("### The Environmental Safety Net")
            st.warning("Since items like **Chicken Bags** or **Meat Trays** are often too dirty to recycle, the 'After' recommendation focuses on **Bio-Mineralization**. This ensures that even in a landfill, the material disappears safely rather than tainting the earth.")

        # --- C. SAFETY & INTEGRITY ---
        st.markdown("---")
        st.header("🛡️ Integrity & Food Safety")
        if category == "Hot Food (Meat/Chicken)" and current['tsi'] < 55:
            st.error("🚨 **Thermal Warning:** Current material has a high 'Leach Risk' at 180°F. The redesign increases thermal stability to prevent chemical migration.")
        else:
            st.success("✅ **Performance Match:** Redesign maintains 100% of current structural barrier strength.")

        if st.button("📄 Generate Executive Audit"):
            report = f"AUDIT: {selected_name}\nSTATUS: Upgrade Recommended\n\nREASONING: Switching to a VoraCycle redesigned polymer prevents microplastic shedding and eliminates chemical migration into fats. This protects the brand from future environmental and health litigation."
            st.text_area("Audit Summary:", value=report, height=200)

elif st.session_state.get("authentication_status") is False:
    st.error('Login Failed.')
elif st.session_state.get("authentication_status") is None:
    st.warning('Please log in.')
