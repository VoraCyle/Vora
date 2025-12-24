import streamlit as st
from openai import OpenAI
import json

# --- 1. ACCESS & DATABASE LOADING ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except:
    st.error("🚨 API Key Missing.")

# Load your Monster Library
try:
    with open('material_library.json', 'r') as f:
        monster_db = json.load(f)
except:
    monster_db = {}

# --- 2. THE AUDIT ENGINE (The 'Before & After' Logic) ---
def run_dna_audit(monster_id, volume):
    if monster_id in monster_db:
        item = monster_db[monster_id]
        savings = item['tax_penalty'] * volume
        
        st.divider()
        st.subheader(f"🧬 DNA Transformation: {item['name']}")
        
        col1, col2 = st.columns(2)
        with col1:
            st.error("**CURRENT DNA (The Monster)**")
            st.write(f"**Material:** {item['description']}")
            st.write(f"**Tax Penalty:** ${item['tax_penalty']} / unit")
        
        with col2:
            st.success("**VORA DNA (The Cure)**")
            st.write(f"**New Material:** {item['vora_fix']['material']}")
            st.write(f"**Annual Savings:** ${savings:,.2f}")

        # The Mixer Prescription
        with st.expander("📝 INDUSTRIAL RECIPE (Send to Mixer)"):
            st.code(item['vora_fix']['recipe'], language="text")
            st.caption("Standardized batch instructions for circular manufacturing.")

# --- 3. UI LAYOUT ---
st.title("🛡️ VoraCycle: Strategic DNA Command")

# Sidebar for Database Selection
st.sidebar.header("Vora Monster Registry")
selected_monster = st.sidebar.selectbox("Select a Material Monster:", ["-- Select --"] + list(monster_db.keys()))

volume = st.sidebar.number_input("Annual Unit Volume:", value=1000000)

if selected_monster != "-- Select --":
    run_dna_audit(selected_monster, volume)
else:
    st.info("Select a 'Monster' from the sidebar to see the DNA Transformation and Mixer Recipe.")