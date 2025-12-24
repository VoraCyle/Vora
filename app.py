import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. THE VORA 100: HARD-TO-IMPROVE ASSET LIBRARY ---
# Categorized for Enterprise Logistics
vora_100 = {
    "Grocery & Food Service": [
        "Rotisserie Chicken Trays (CPET/APET)", "Multi-Layer Snack Pouches (Mylar)", 
        "Black Plastic Meat Trays", "Waxed Produce Boxes", "K-Cups / Coffee Pods",
        "BOPP Bread Bags", "Frozen Food Poly-bags", "Polystyrene Egg Cartons",
        "Aseptic Juice Cartons (Tetra Pak)", "Thermal Receipt Paper (BPA/BPS)"
    ],
    "Health, Beauty & Pharmacy": [
        "Multi-Material Blister Packs", "Mixed-Plastic Squeeze Tubes", 
        "Pump-Action Dispenser Heads", "Non-Recyclable Makeup Compacts",
        "Vitamins/Supplement Bottles (Colored PET)", "Aerosol Sunscreen Cans",
        "Sheet Mask Sachet Laminates", "Plastic Cotton Swab Stems"
    ],
    "Warehouse & Logistics": [
        "LLDPE Stretch Wrap (Dirty/Stretched)", "Nylon Strapping Bands", 
        "EPS (Styrofoam) Packing Peanuts", "Chemically Treated Timber Pallets",
        "Bubble Mailers (Mixed Paper/Plastic)", "Heavy-Duty Vinyl Tarps",
        "Adhesive Shipping Labels (Non-Recyclable)", "Polypropylene Shipping Totes"
    ],
    "Hardlines & Electronics": [
        "Lithium-Ion Battery Slurry", "Integrated LED Light Strips",
        "Mixed-Metal Power Cables", "Non-Removable Battery Tools",
        "PVC Clamshell Packaging", "Carbon Fiber Composite Scrap",
        "Small Appliance Plastic Housings", "Remote Control Assemblies"
    ],
    "Textiles & Apparel": [
        "Spandex/Cotton Blend Clothing", "Synthetic Microfiber Fleece",
        "Multi-Material Athletic Shoes", "Plastic Hanger Waste",
        "Polybag Clothing Protection", "Treated Waterproof Outerwear"
    ]
}

# Flatten for the dropdown
dropdown_items = ["-- Select a 'Hard 100' Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are the VoraCycle Lead Arbiter. You provide technical engineering blueprints and supplier grading for circular retail systems. You focus on 'Start-at-the-Beginning' DNA optimization."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: The Hard 100", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.markdown("### Engineering for 'Universal Sustainability'—Safe in Any Waste Stream.")

col1, col2 = st.columns(2)
with col1:
    st.markdown("##### 📦 1. High-Volume 'System Spoilers'")
    dropdown_choice = st.selectbox("Select a problematic retail asset:", dropdown_items)

with col2:
    st.markdown("##### 🔍 2. Custom SKU / Concept Audit")
    search_query = st.text_input("Type in any specific product or material:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a 'Hard 100' Strategic Asset --" else None)

if final_query:
    st.divider()
    with st.spinner(f"Simulating Forensic Path for {final_query}..."):
        
        master_prompt = f"""
        Conduct a Forensic DNA Audit for: {final_query}.
        
        ### 📊 DUAL-PATH FORENSIC SCORECARD
        Compare 'Status Quo' vs. 'Vora Resilient Design'.
        Include a Table: Primary Outcome, Sustainability Rating, Capital Retention %, and Environmental Toxin Risk.

        ### 🧬 DNA ENGINEERING PROTOCOL: THE "WHAT & HOW"
        Provide a technical blueprint for the design team:
        1. **WHAT TO REMOVE:** Identify specific toxins or non-separable materials.
        2. **WHAT TO REPLACE WITH:** Suggest mono-material or bio-mineral alternatives.
        3. **HOW TO ASSEMBLE:** Describe the 'Design for Disassembly' or 'Molecular Purity' process needed so the item is safe for the planet even if it enters a landfill. (300+ Words).

        ### 🎯 THE BEST PATH & HANDLING
        Clearly identify the #1 most profitable path. Provide step-by-step handling and backhauling logistics for a company like Costco.

        ### 📋 VORA SUPPLIER SCORECARD
        Evaluate a theoretical supplier providing this item based on the new DNA.
        Format as a Table: Material Purity, Chemical Transparency, Circular Recovery Readiness, and Overall Vora Grade (A-F).

        ---

        ### 🏁 FINAL EXECUTIVE SUMMARY: THE SUSTAINABILITY LEADERSHIP VERDICT
        Why is this redesign the only logical choice for the company's 10-year survival? How does it eliminate consumer liability? (200 Words).
        """

        st.markdown(generate_vora_analysis(master_prompt))

# --- 5. FOOTER ---
st.sidebar.info(f"VoraCycle v4.9.0 | {sum(len(v) for v in vora_100.values())} System Spoilers Cataloged") 