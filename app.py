import streamlit as st
from openai import OpenAI
import json

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. MATERIAL LIBRARY LOADING ---
# This pulls the "Monster Book" into the app's memory
try:
    # Ensure your file is named exactly 'material_library.json' in your folder
    with open('material_library.json', 'r') as f:
        monster_db = json.load(f)
except FileNotFoundError:
    st.error("🚨 'material_library.json' not found. Please ensure the file is in your project folder.")
    monster_db = {}

# --- 3. VORA AUDIT LOGIC ---
def run_audit(monster_id, volume):
    """
    Analyzes a material monster and returns the Vora DNA recipe for the mixer.
    """
    if monster_id in monster_db:
        item = monster_db[monster_id]
        savings = item['tax_penalty'] * volume
        
        # Display results in the app
        st.success(f"✅ Vora Fix: {item['vora_fix']['material']}")
        st.metric("Potential Annual Savings", f"${savings:,.2f}")
        
        # This is the 'Prescription' to be sent to the Mixer/Blender
        return item['vora_fix']['recipe']
    return None

# --- 4. STREAMLIT UI ---
st.title("Vora DNA: Forensic Material Audit")
st.markdown("Bridge the gap between **Retail Liability** and **Manufacturing Solutions**.")

# Sidebar for Library Selection
st.sidebar.header("Vora Monster Registry")
monster_options = list(monster_db.keys())
selected_id = st.sidebar.selectbox("Select a Material Monster:", ["-- Select --"] + monster_options)

# User Inputs
col1, col2 = st.columns(2)
with col1:
    volume = st.number_input("Annual Unit Volume:", min_value=0, value=1000000, step=100000)
with col2:
    st.info("Calculates 'Material Debt' based on regional tax penalties.")

# Run Audit Display
if selected_id != "-- Select --":
    st.divider()
    # Execute the audit and capture the recipe for the mixer
    recipe = run_audit(selected_id, volume)
    
    if recipe:
        st.subheader("🧬 Vora DNA Prescription")
        st.write("Send the following technical specs to your third-party pellet mixer:")
        
        with st.expander("📝 View Mixer Recipe (Technical Data)"):
            st.code(recipe, language="text")
            st.caption("Batch-ready instruction for industrial extrusion lines.")

# Custom AI Audit
st.divider()
search_query = st.text_input("Custom SKU/Material Audit (AI Analysis):")

if search_query:
    with st.spinner(f"Analyzing {search_query}..."):
        # Existing OpenAI prompting logic goes here
        st.write("AI Forensic Analysis complete.")