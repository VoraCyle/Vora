import streamlit as st
from openai import OpenAI
import json

# --- 1. ACCESS & DATABASE ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except:
    st.error("🚨 API Key Missing.")

# Load the Monster Library (Memory)
try:
    with open('material_library.json', 'r') as f:
        monster_db = json.load(f)
except:
    monster_db = {}

# --- 2. THE VORA 100 STRATEGIC REGISTRY ---
vora_100 = {
    "🥩 POULTRY & FRESH MEATS": ["MAP Poultry Trays", "Absorbent Poultry Pads", "Vacuum Wraps", "Black Meat Trays"],
    "🧻 PAPER & HYGIENE WRAPS": ["Kirkland Bath Tissue Case-Wrap", "Paper Towel Individual Rolls"],
    "❄️ FROZEN & REFRIGERATED": ["Aqueous Frozen Bags", "Multi-Layer Meal Pouches"],
    "🚩 HIGH-RISK LIABILITIES": ["PVC Clamshells", "PFAS Wrappers", "LLDPE Stretch Wrap"]
}

# --- 3. THE FORENSIC ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are the VoraCycle Chief Sustainability Officer. Focus on Pre-Emptive Circularity and Industrial DNA Reconstruction."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.4 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Executive Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.markdown("### The Pre-Emptive Circularity Report: Future-Proofing the Enterprise.")

# Build Dropdown List
dropdown_items = ["-- Select a Strategic Asset --"]
for category, items in vora_100.items():
    dropdown_items.extend(items)

col1, col2 = st.columns(2)
with col1:
    choice = st.selectbox("Select Asset for Executive Audit:", dropdown_items)
with col2:
    search = st.text_input("Search Custom SKU (e.g. 'Kirkland Salmon'):")

volume = st.sidebar.number_input("Annual Unit Volume:", value=1000000)
final_query = search if search else (choice if choice != "-- Select a Strategic Asset --" else None)

if final_query:
    st.divider()
    
    # Check if we have HARD DATA for this in the JSON
    if final_query in monster_db:
        item = monster_db[final_query]
        savings = item['tax_penalty'] * volume
        st.success(f"✅ HARD DATA MATCH: {item['name']}")
        st.metric("Potential Annual Savings", f"${savings:,.2f}")
        
        with st.expander("📝 INDUSTRIAL DNA RECIPE (Send to Mixer)"):
            st.code(item['vora_fix']['recipe'], language="text")
            st.caption("Technical instructions for standardized industrial extrusion.")

    # Generate the Executive Impact Report
    with st.spinner(f"Generating Executive Impact Report for {final_query}..."):
        master_prompt = f"Generate a Pre-Emptive Impact Report for: {final_query}. Include: 1. Executive Dashboard (BAU vs Vora), 2. Ecological Hedge, 3. Financial Fortress, 4. DNA Failsafe."
        st.markdown(generate_vora_analysis(master_prompt))