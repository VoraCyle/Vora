import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. UPDATED COSTCO-SPECIFIC ASSET REGISTRY ---
vora_100 = {
    "🍗 THE ROTISSERIE CATEGORY": [
        "Multi-Layer Chicken Bags (High-Heat Composite)", 
        "Absorbent Chicken Pad (Poly-Fiber)", 
        "Rigid Plastic Chicken Trays (Legacy SKU)",
        "PP Chicken Bag Handles"
    ],
    "🧻 PAPER & HYGIENE WRAPS": [
        "Kirkland Bath Tissue Case-Wrap (LDPE)", 
        "Paper Towel Individual Rolls (LDPE Thin-Film)", 
        "Facial Tissue Box Windows (Film/Fiber)",
        "Napkin Bundle Overwrap"
    ],
    "🥩 POULTRY & FRESH MEATS": [
        "MAP Poultry Trays", "Absorbent Poultry Pads", "Vacuum Wraps", "Black Meat Trays"
    ],
    "❄️ FROZEN & REFRIGERATED": [
        "Aqueous Frozen Bags", "Multi-Layer Meal Pouches", "Mono-PE Trays", "Poly-Ice Cream Cartons"
    ],
    "🚩 HIGH-RISK LIABILITIES": [
        "PVC Clamshells", "PFAS Wrappers", "Lithium Battery Packs", "LLDPE Warehouse Stretch Wrap"
    ]
}

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": """You are the VoraCycle Lead Forensic Engineer. 
                Focus on 'Path-Agnostic' DNA. For high-volume items like chicken bags and tissue wraps, 
                explain how to move from multi-material composites to mono-material or bio-mineralized 
                DNA that is 100% food-safe and landfill-safe."""},
                {"role": "user", "content": prompt}
            ],
            temperature=0.4 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: High-Volume Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.markdown("### Engineering the 'Hard 100': From Chicken Bags to Case-Wraps.")

dropdown_items = ["-- Select a Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

col1, col2 = st.columns(2)
with col1:
    dropdown_choice = st.selectbox("Select high-volume Warehouse SKU:", dropdown_items)
with col2:
    search_query = st.text_input("Custom SKU/Material Audit:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a Strategic Asset --" else None)

if final_query:
    st.divider()
    with st.spinner(f"Analyzing {final_query} for Pre-Emptive ROI..."):
        
        master_prompt = f"""
        Execute a Forensic DNA Audit for: {final_query}.

        ### 📊 DUAL-PATH FORENSIC SCORECARD
        Compare 'Status Quo' (Current Design) vs. 'Vora Resilient Design'.

        ### 🧬 THE DNA RECONSTRUCTION (Universal Safety)
        1. **PATH A (Waste/Trash):** If this {final_query} ends up in a trash bin, how does the new DNA ensure it doesn't leak toxins and mineralizes safely?
        2. **PATH B (Recycle/Loop):** How do we simplify the film (e.g. tissue wrap) to be a high-value mono-PE that can be recycled back into our own warehouse bags?

        ### 🍗 FOOD SAFETY & PURITY (If Applicable)
        Confirm that high-heat resistance for items like chicken bags is maintained without using PFAS or toxic glues.

        ### 💰 COST-BENEFIT ESTIMATOR (The Costco Logic)
        - **EPR Fee Reduction:** How much do we save in plastic taxes by switching from composite to mono-material?
        - **Logistics Gain:** How does the new DNA reduce weight or shipping volume?

        ### 🏁 EXECUTIVE VERDICT
        Why is this rebuild essential for Kirkland Signature brand trust?
        """

        st.markdown(generate_vora_analysis(master_prompt))