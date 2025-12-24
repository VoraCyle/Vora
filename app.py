
import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. EXPANDED STRATEGIC ASSET REGISTRY ---
vora_100 = {
    "🥩 POULTRY & FRESH MEATS": [
        "MAP (Modified Atmosphere) Poultry Trays", "Absorbent Poultry Pads (Composite)", 
        "Vacuum-Sealed Poly-Vinyl Wraps", "Black Plastic Meat Trays", 
        "Soaker Pads (Chemical-Bound Fiber)"
    ],
    "🥬 FRESH PRODUCE & GOODS": [
        "Cellulose Berry Clamshells", "Breathable Mixed-Poly Produce Bags", 
        "Waxed Fruit Corrugated Boxes", "Plastic Mesh Citrus Bags",
        "Herbal Sachet Laminates"
    ],
    "❄️ FROZEN & REFRIGERATED": [
        "Aqueous-Coated Frozen Vegetable Bags", "Multi-Layer Frozen Meal Pouches", 
        "High-Barrier Mono-PE Trays", "Poly-Lined Ice Cream Cartons",
        "Chilled Liquid Cartons (Aluminum-Poly)"
    ],
    "📦 DRY GOODS & PANTRY": [
        "Metallized PP Snack Liners", "Composite Cardboard Canisters", 
        "BOPP Cereal Liners", "Multi-Wall Pet Food Bags",
        "Plastic-Spout Pour Bags"
    ],
    "🚩 HIGH-RISK LIABILITIES": [
        "PVC Clamshell Packs", "PFAS-Coated Fast Food Wrappers", 
        "Lithium-Ion Tool Battery Packs", "Industrial LLDPE Stretch Wrap",
        "Thermal BPA Receipt Paper"
    ]
}

# Flatten for dropdown
dropdown_items = ["-- Select a Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": """You are the VoraCycle Lead Forensic Engineer. 
                Focus: Path-Agnostic DNA Engineering. 
                Critical Requirement: All DNA transformations for food-contact items MUST 
                be certified Bio-Inert and Food-Safe. Ensure the reconstruction eliminates 
                the risk of toxin migration or chemical leaching into food products."""},
                {"role": "user", "content": prompt}
            ],
            temperature=0.4 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Food-Safe DNA", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.markdown("#### Engineering for 'Universal Sustainability' & Food-Contact Safety.")

col1, col2 = st.columns(2)
with col1:
    st.markdown("##### 📦 1. Strategic Asset Selection")
    dropdown_choice = st.selectbox("Choose a known retail vulnerability:", dropdown_items)

with col2:
    st.markdown("##### 🔍 2. Custom SKU / Material Audit")
    search_query = st.text_input("Search custom food packaging or SKU:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a Strategic Asset --" else None)

if final_query:
    st.divider()
    with st.spinner(f"Engineering Food-Safe DNA for {final_query}..."):
        
        master_prompt = f"""
        Conduct a Path-Agnostic Forensic Audit for: {final_query}.

        ### 📊 DUAL-PATH FORENSIC SCORECARD
        Compare 'Status Quo' vs. 'Vora Resilient Design'.

        ### 🍎 FOOD SAFETY & DNA PURITY PROTOCOL
        Explain how this new DNA protects the food:
        - **Non-Toxic Barrier:** How do we eliminate PFAS, toxic glues, and chemical leaching?
        - **Bio-Inert Integrity:** Explain why the new material (e.g., Aqueous coatings or PHAs) is safer for food contact than the status quo.
        - **Temperature Stability:** Confirm safety for frozen, chilled, or ambient storage.

        ### 🧬 THE DNA FAILSAFE (Path-Agnostic)
        - **PATH A (Waste/Landfill):** How does the DNA ensure safe mineralization?
        - **PATH B (Recycle/Technical):** How do we simplify the DNA for 100% recovery?

        ### 💰 FINANCIAL ROI & LIABILITY
        Estimate EPR fine savings and 'Brand Trust' value for a company like Costco.

        ### 📋 VORA SUPPLIER SCORECARD
        Evaluate the manufacturer. Table: Material Purity, Chemical Transparency, Circular Readiness, Vora Grade.
        """

        st.markdown(generate_vora_analysis(master_prompt))
else:
    st.info("👆 Select a food packaging category or type a custom material to begin.")

# --- 5. FOOTER ---
st.sidebar.info(f"VoraCycle v5.5.0 | Total Food-Safe DNA Profiles: {sum(len(v) for v in vora_100.values())}")