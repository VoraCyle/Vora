import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. COSTCO-SCALE "HARD" ASSETS ---
problematic_items = [
    "-- Select a High-Volume Asset --",
    "Multi-Layer Pet Food Bags (Mylar/Poly)",
    "Black Plastic Food Trays (Carbon Black)",
    "LLDPE Pallet Shrink Wrap",
    "Expanded Polystyrene (EPS) Coolers",
    "Lithium-Ion Tool Batteries (Returns)",
    "Wax-Coated Produce Corrugated Boxes",
    "PVC Blister Packaging",
    "Composite Flooring Returns",
    "Textile Waste (Kirkland Clothing)",
    "Industrial Pallet Strapping (Nylon/PP)"
]

# --- 3. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are the VoraCycle Lead Arbiter. You provide forensic comparative audits for enterprise retailers like Costco. You focus on design-phase DNA optimization and logistics."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Forensic Command", layout="wide")
st.title("🛡️ VoraCycle: Strategic Command Center")
st.markdown("#### Comparative Forensic Audit: Status Quo vs. Vora Optimized")

col1, col2 = st.columns(2)
with col1:
    dropdown_choice = st.selectbox("Select problematic Warehouse Stream:", problematic_items)
with col2:
    search_query = st.text_input("Custom SKU/Material Audit:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a High-Volume Asset --" else None)

if final_query:
    st.divider()
    with st.spinner("Executing Dual-Path Simulation..."):
        
        master_prompt = f"""
        Audit the following for an Enterprise Retailer: {final_query}.

        ### 📊 DUAL-PATH FORENSIC SCORECARD
        Compare the 'Status Quo' (Current Disposal/Recycling) vs. 'Vora Optimized' (Highest Circularity).
        Output this as a Markdown Table with the following rows:
        - Primary Outcome (e.g. Landfill vs Refurbish)
        - Sustainability Rating (0-10)
        - Capital Retention (%)
        - Environmental Impact Score (0-10)
        - Risk/Liability Level (High/Low)

        ---

        ### ⚖️ THE WINNING VERDICT: BEST STRATEGIC PATH
        Clearly state the chosen Vora Path. Provide 200 words of forensic logic explaining why this is the best for the ENVIRONMENT and the COMPANY. 

        ### 🚚 HANDLING & LOGISTICS (SOP)
        1. **Warehouse Handling:** Precise sorting/prepping instructions.
        2. **Logistics:** Transport path details and 'Backhauling' opportunities.
        3. **Facility:** Specify the target Circular Hub or Processing center.

        ### 🧬 DNA MODIFICATION (THE DESIGN FIX)
        How must we re-engineer the PRODUCT DNA right now? Detail how to make it 'Safe by Design' so it acts as a non-toxic biological nutrient if it ever reaches end-of-life. (250+ Words).

        ---

        ### 🏁 FINAL EXECUTIVE SUMMARY: THE VORA ADVANTAGE
        Provide a 200-word concluding argument. 
        1. WHY: Why is the Vora Optimized path the most beneficial for the long-term health of the business? 
        2. HOW: How does this specific strategy keep sustainability high while eliminating waste before it starts?
        3. SYNERGY: Explain the 'Zero-Leakage' promise—where profit and planet align perfectly.
        """

        st.markdown(generate_vora_analysis(master_prompt))

# --- 5. FOOTER ---
st.sidebar.info("VoraCycle Enterprise Logic: v4.4.0")