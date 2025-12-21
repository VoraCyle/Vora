import streamlit as st
import google.generativeai as genai
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
import pandas as pd

# --- 1. SECURE CONFIGURATION ---
if "GEMINI_API_KEY" in st.secrets:
    genai.configure(api_key=st.secrets["GEMINI_API_KEY"])
    model = genai.GenerativeModel('gemini-1.5-flash')
else:
    st.error("🔑 API Key Missing. Please add GEMINI_API_KEY to Streamlit Secrets.")
    st.stop()

# --- 2. GLOBAL MATERIAL INVENTORY ---
product_inventory = {
    "Search or select an item...": "",
    "Meat Wrap (PVC)": "C=CCl",
    "Water Bottle (PET)": "CC1=CC=C(C=C1)C(=O)OCCO",
    "Chip Bag (BOPP/Multi-layer)": "CCCCCCCCCC.C=CC#N",
    "Deli Container (PP)": "CC(C)CC(C)C",
    "Frozen Food Bag (LDPE)": "CCCCCCCCCCCC",
    "Pharmacy Bottle (PC)": "CC(C)(C1=CC=C(OC(=O)OC2=CC=C(C(C)(C)C)C=C2)C=C1)C",
    "Coffee Cup Liner (PE)": "CCCCCCCC",
}

# --- 3. THE ANALYTICAL DECISION ENGINE ---
def run_synthesis_audit(item_name, smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        toxic = any(a.GetSymbol() in ['Cl', 'F', 'Br', 'I'] for a in mol.GetAtoms())
        
        # Path 1: Recycle Logic
        r_rating = 92 if ("PET" in item_name or "Bottle" in item_name) else 38
        
        # Path 2: Landfill Logic
        l_rating = 15 if toxic else 42
        
        # FINAL INTEGRATED RATING (Weighted average of current reality)
        # Most items are more likely to hit landfill, so we weigh that risk
        final_before_score = round((r_rating * 0.3) + (l_rating * 0.7), 1)
        
        return r_rating, l_rating, final_before_score, toxic
    except:
        return None

# --- 4. THE APEX INTERFACE ---
st.set_page_config(page_title="VoraCycle Apex OS", layout="wide")
st.title("🔮 Wraith VoraCycle: Synthesis Apex OS")
st.markdown("### *Integrated Endgame Audit: From Liability to Asset*")

search_query = st.selectbox("🧬 Select Product for Forensic Audit:", list(product_inventory.keys()))

if search_query and search_query != "Search or select an item...":
    active_smiles = product_inventory[search_query]
    audit = run_synthesis_audit(search_query, active_smiles)
    
    if audit:
        r_score, l_score, before_total, is_toxic = audit
        after_total = 98.4  # VoraCycle Standard for Grade A Mineralization
        
        # --- THE START-LINE RATINGS ---
        st.divider()
        st.header("📊 Integrated Forensic Scores")
        col_metrics = st.columns(3)
        col_metrics[0].metric("Legacy Score (Before)", f"{before_total}%", delta="CRITICAL RISK" if before_total < 50 else "STABLE")
        col_metrics[1].metric("VoraCycle Score (After)", f"{after_total}%", delta="GRADE A TARGET")
        col_metrics[2].metric("Best Pathway", "VoraCycle Mineralization", delta="100% Circular")

        # --- THE PATHWAY COMPARISON ---
        st.divider()
        path_a, path_b = st.tabs(["🚀 PATH A: RECYCLING LIMITS", "🌋 PATH B: LANDFILL IMPACT"])
        
        with path_a:
            st.subheader(f"Mechanical Recycling Potential: {r_score}%")
            st.write("""
            **Summary:** Most legacy materials lose structural integrity after 2-3 recycling cycles. 
            For food-grade items (frozen/fresh), 'Downcycling' is common, meaning the product 
            is turned into low-value items like park benches rather than new packaging.
            """)
            

        with path_b:
            st.subheader(f"Landfill Environmental Impact: {l_score}%")
            st.write(f"""
            **Summary:** Currently, this material is a long-term liability. In landfill conditions, 
            it stays mummified or leaches toxins. For fresh food packaging, this leads to 
            microplastic shedding that enters the food chain.
            """)

        # --- THE DEEP SUMMARY (AI REASONING) ---
        st.divider()
        st.header("⚖️ Executive Transformation Summary")
        
        prompt = (
            f"Provide a deep forensic summary for the product '{search_query}'. "
            f"1. Why are the current paths (Recycle/Landfill) failing? "
            f"2. Detail the exact molecular changes needed (Introduction of Metabolic Handles). "
            f"3. Explain how these changes improve the product without affecting fresh, frozen, or dry food. "
            f"4. Justify why VoraCycle Mineralization is the superior path for a business like Costco."
        )
        
        with st.spinner("Synthesizing forensic data..."):
            try:
                response = model.generate_content(prompt)
                st.info(response.text)
            except Exception as e:
                st.error(f"📡 API Connection Error: {str(e)}")

        # --- THE "HOW & WHY" VISUALIZATION ---
        st.divider()
        st.subheader("🛠️ The VoraCycle Surgery: How & Why")
        col_how, col_why = st.columns(2)
        
        with col_how:
            st.markdown("**How the changes are made:**")
            st.write("""
            We modify the polymer backbone by inserting 'Latent Scission Bridges.' 
            These bridges are chemically inactive while holding frozen or fresh food. 
            They only break when the 'Enzymatic Key' (found in soil bacteria) is applied.
            """)
            
            
        with col_why:
            st.markdown("**Why this is the superior path:**")
            st.write("""
            Traditional recycling is inconsistent. Landfills are expensive. 
            By making the material 'Soil-Safe' at the start-line, we ensure 
            that no matter where the package ends up, it returns to the earth 
            as water and CO2. This protects the business from future plastic taxes.
            """)
            

    else:
        st.warning("Forensic analysis failed. Please check the material library.")
