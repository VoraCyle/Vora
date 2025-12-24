import streamlit as st
from openai import OpenAI
import json

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except:
    st.error("🚨 API Key Missing.")

# Load the Monster Library
try:
    with open('material_library.json', 'r') as f:
        monster_db = json.load(f)
except:
    monster_db = {}

# --- 2. THE VORA 100 STRATEGIC REGISTRY ---
vora_100 = {
    "🥩 POULTRY & MEATS": ["Multi-Layer Chicken Bags", "MAP Poultry Trays", "Absorbent Poultry Pads", "Black Meat Trays"],
    "🧻 HYGIENE WRAPS": ["Kirkland Bath Tissue Case-Wrap", "Paper Towel Overwrap"],
    "🚩 HIGH-RISK": ["PVC Clamshells", "PFAS Wrappers", "LLDPE Stretch Wrap"]
}

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[{"role": "system", "content": "You are the VoraCycle CSO. Focus on DNA Engineering."}],
            temperature=0.4 
        )
        return response.choices[0].message.content
    except:
        return "Analysis Error."

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: DNA Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")

# Sidebar for Logistics
st.sidebar.header("📦 Production Logistics")
volume = st.sidebar.number_input("Annual Unit Volume:", value=1000000)
unit_weight_grams = st.sidebar.slider("Avg. Unit Weight (grams):", 5, 100, 25)

dropdown_items = ["-- Select Asset --"]
for cat, items in vora_100.items():
    dropdown_items.extend(items)

col1, col2 = st.columns(2)
with col1:
    choice = st.selectbox("Select Asset for Executive Audit:", dropdown_items)
with col2:
    search = st.text_input("Search Custom SKU:")

final_query = search if search else (choice if choice != "-- Select Asset --" else None)

if final_query:
    st.divider()
    
    # --- 5. THE NEW DNA CALCULATOR (WORKS FOR ALL) ---
    st.subheader(f"🧪 New DNA Recipe: {final_query}")
    
    # Calculate Total Tonnage
    total_weight_kg = (volume * unit_weight_grams) / 1000
    total_tons = total_weight_kg / 1000
    
    st.metric("Total Material Volume Required", f"{total_tons:,.2f} Metric Tons")

    # Get recipe from DB or use a Default AI Template
    if final_query in monster_db:
        recipe_raw = monster_db[final_query]['vora_fix']['recipe']
    else:
        # Default Vora DNA Template for unknown items
        recipe_raw = "92% Mono-Polymer Base, 5% Vora-C1 Catalyst, 3% Bio-Mineral Tracer"

    recipe_parts = recipe_raw.split(',')
    cols = st.columns(len(recipe_parts))
    
    for idx, part in enumerate(recipe_parts):
        try:
            # Extract number before %
            percent_str = part.split('%')[0].strip()
            percent = float(''.join(filter(lambda x: x.isdigit() or x=='.', percent_str))) / 100
            part_weight = total_tons * percent
            with cols[idx]:
                st.info(f"**DNA Element {idx+1}**\n\n{part.strip()}\n\n**Order:** {part_weight:,.2f} Tons")
        except:
            with cols[idx]:
                st.info(f"**DNA Element {idx+1}**\n\n{part.strip()}")
    
    st.caption("ℹ️ *This manifest represents the required 'New DNA' components for the industrial mixer.*")
    st.divider()

    # 6. EXECUTIVE REPORT
    with st.spinner(f"Analyzing {final_query}..."):
        master_prompt = f"Generate a Pre-Emptive Impact Report for: {final_query}."
        st.markdown(generate_vora_analysis(master_prompt))