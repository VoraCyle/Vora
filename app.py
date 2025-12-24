import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. THE VORA 100 STRATEGIC REGISTRY ---
vora_100 = {
    "🥩 POULTRY & FRESH MEATS": ["MAP Poultry Trays", "Absorbent Poultry Pads", "Vacuum Wraps", "Black Meat Trays"],
    "🥬 FRESH PRODUCE & GOODS": ["Cellulose Berry Clamshells", "Bio-Produce Bags", "Waxed Boxes", "Mesh Citrus Bags"],
    "❄️ FROZEN & REFRIGERATED": ["Aqueous Frozen Bags", "Multi-Layer Meal Pouches", "Mono-PE Trays", "Poly-Ice Cream Cartons"],
    "📦 DRY GOODS & PANTRY": ["Metallized Snack Liners", "Composite Canisters", "BOPP Cereal Liners", "Multi-Wall Pet Food Bags"],
    "🚩 HIGH-RISK LIABILITIES": ["PVC Clamshells", "PFAS Wrappers", "Lithium Battery Packs", "LLDPE Stretch Wrap", "BPA Receipts"]
}

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": """You are the VoraCycle Lead Forensic Engineer. 
                Your focus is 'Pre-Emptive Circularity'. You prove that fixing the end-game at the start 
                is the ultimate financial and ecological hedge. 
                Explain how fixing the DNA creates a 'Universal Success' outcome where the company 
                can never lose, regardless of what the consumer does with the item."""},
                {"role": "user", "content": prompt}
            ],
            temperature=0.4 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Pre-Emptive Circularity", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.markdown("### Deleting Waste Before It Starts: The Pre-Emptive Circularity Model.")

dropdown_items = ["-- Select a Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

col1, col2 = st.columns(2)
with col1:
    dropdown_choice = st.selectbox("Select a known retail vulnerability:", dropdown_items)
with col2:
    search_query = st.text_input("Search custom food packaging or SKU:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a Strategic Asset --" else None)

if final_query:
    st.divider()
    with st.spinner(f"Simulating Pre-Emptive ROI for {final_query}..."):
        
        master_prompt = f"""
        Execute a Pre-Emptive DNA Audit for: {final_query}.

        ### 🌍 THE ECOLOGICAL EVOLUTION: DELETING THE LANDFILL
        - **DNA Fix vs. Status Quo:** How does fixing the DNA at the start prevent microplastics and forever chemicals (PFAS)?
        - **The 'Universal Safety' Outcome:** Explain how the product becomes 'Earth-Native'. If it ends up in the ocean or soil, how does this new DNA allow it to mineralize safely?
        - **Carbon Reduction:** How does a simplified DNA reduce the energy required for manufacturing and transport?

        ### 💰 THE FINANCIAL FORTRESS: PROFIT THROUGH PURITY
        - **EPR Penalty Immunity:** How does this design make the company immune to 'Plastic Taxes' and 'Waste Fines' (California SB 54, etc.)?
        - **Reverse Logistics Profit:** How does a pure DNA allow the company to buy back its own waste as raw material, lowering future costs?
        - **Supply Chain Resilience:** How does removing 'Material Monsters' simplify sourcing and protect against price spikes?

        ### 🧬 THE DNA FAILSAFE (Path-Agnostic)
        - **PATH A (Waste):** Safe mineralization protocol.
        - **PATH B (Recycle):** 100% technical recovery protocol.
        
        ### 🍎 FOOD SAFETY & DNA PURITY
        - Confirm non-toxic migration and bio-inert integrity.

        ### 🏁 EXECUTIVE VERDICT: THE FUTURE-PROOF COMPANY
        Explain how this makes the company the #1 sustainability leader by removing the 'burden of choice' from the consumer.
        """

        st.markdown(generate_vora_analysis(master_prompt))
else:
    st.info("👆 Start by selecting an asset to see the Pre-Emptive Circularity model in action.")

# --- 5. FOOTER ---
st.sidebar.info(f"VoraCycle v5.6.0 | Pre-Emptive Circularity Engine")