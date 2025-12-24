import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])

# --- 2. STRATEGIC ASSET CATEGORIES ---
vora_100 = {
    "🚩 TOP LIABILITIES (High EPR Fines)": [
        "Multi-Layer Mylar Snack Pouches", "Black Plastic Rotisserie Trays", 
        "PVC Clamshell Electronics Packs", "PFAS-Coated Fast Food Wraps",
        "Aseptic Juice Cartons (Tetra Pak)", "BOPP Non-Recyclable Pet Food Bags"
    ],
    "⚡ QUICKEST FIXES (Immediate ROI)": [
        "LLDPE Warehouse Stretch Wrap", "BPA-Free Thermal Receipts", 
        "Mono-material Polyethylene Mailers", "Uncoated Corrugated Cardboard",
        "Natural Fiber Pallet Strapping", "Detectable Pigment Black Plastics"
    ],
    "🧬 LONG-TERM DNA SHIFTS (High Impact)": [
        "Lithium-Ion Battery Slurry", "Mixed-Textile Kirkland Apparel",
        "Integrated LED Appliance Panels", "Multi-Material Athletic Footwear",
        "EPS (Styrofoam) Cold-Chain Coolers", "Composite Construction Returns"
    ]
}

# --- 3. THE "REAL-WORLD" ENGINE ---
def generate_vora_analysis(prompt):
    # We use 'gpt-4o' for maximum reasoning power to avoid "template" answers
    response = client.chat.completions.create(
        model="gpt-4o", 
        messages=[
            {"role": "system", "content": """You are the VoraCycle Lead Forensic Engineer. 
            DO NOT give generic 'mono-material' advice unless it is the only viable path. 
            For every item, analyze its specific chemical and mechanical makeup. 
            Provide deep-tech solutions: Molecular recycling for polymers, Pyrolysis for composites, 
            Hydrometallurgy for batteries, and Bio-polymers for food-contact items. 
            Your goal is the 'End-Game'—the absolute best path for the planet and P&L."""},
            {"role": "user", "content": prompt}
        ],
        temperature=0.4 # Lower temperature ensures more factual, engineering-heavy results
    )
    return response.choices[0].message.content

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Forensic Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic Liability Command")

# (Dropdown and Search logic remains the same...)
dropdown_items = ["-- Select a Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

col1, col2 = st.columns(2)
with col1:
    dropdown_choice = st.selectbox("Select Asset Category:", dropdown_items)
with col2:
    search_query = st.text_input("Search Custom SKU / DNA Audit:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a Strategic Asset --" else None)

if final_query:
    st.divider()
    with st.spinner(f"Performing Deep Forensic Audit on {final_query}..."):
        
        master_prompt = f"""
        Execute a 2025 Forensic Audit for: {final_query}.
        
        ### 📊 DUAL-PATH COMPARISON
        Show a Markdown Table comparing 'Status Quo' vs. 'Vora Optimized'.
        
        ### ⚙️ THE RECONSTRUCTION BLUEPRINT
        What is the specific mechanical or chemical reconstruction needed? 
        If it's a battery, talk about cathode recovery. If it's a textile, talk about fiber-to-fiber recycling. 
        Provide the 'How'—the exact facility type and process.

        ### 💰 FINANCIAL ROI & LIABILITY
        Estimate savings in EPR fines (e.g., California SB 54 or EU Packaging Directives). 
        How does this redesign lower the 'Total Cost of Ownership'?

        ### 🧬 DNA TRANSFORMATION (The 'Before' Fix)
        How do we change the DNA at the start? List specific chemical swaps (e.g., replace PVC with PE, or remove PFAS for Aqueous coatings).

        ### 📋 SUPPLIER SCORECARD
        Grade a supplier based on their ability to implement this SPECIFIC engineering fix.
        """

        st.markdown(generate_vora_analysis(master_prompt))