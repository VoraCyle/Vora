import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. THE ASSET REGISTRY ---
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
                {"role": "system", "content": "You are the VoraCycle Lead Arbiter. You specialize in 'Resilient Circularity'—engineering products so they remain sustainable regardless of how the consumer disposes of them."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7 
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 4. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle: Resilient Design", layout="wide")
st.title("🛡️ VoraCycle: Strategic DNA Command Center")
st.markdown("### Engineering for 'Universal Sustainability'—Safe in Any Waste Stream.")

col1, col2 = st.columns(2)
with col1:
    dropdown_choice = st.selectbox("Select problematic Warehouse Stream:", problematic_items)
with col2:
    search_query = st.text_input("Custom SKU/Material Audit:")

final_query = search_query if search_query else (dropdown_choice if dropdown_choice != "-- Select a High-Volume Asset --" else None)

if final_query:
    st.divider()
    with st.spinner("Simulating Universal Sustainability Outcomes..."):
        
        master_prompt = f"""
        Audit the following for an Enterprise Retailer: {final_query}.

        ### 📊 DUAL-PATH FORENSIC SCORECARD
        Compare 'Status Quo Design' vs. 'Vora Resilient Design'.
        - Sustainability Rating (0-10)
        - Outcome if Landfilled (Toxin vs Nutrient)
        - Outcome if Recycled (Downcycled vs Upcycled)
        - Financial Risk (High vs Zero)

        ---

        ### 🧬 DNA TRANSFORMATION: THE "START-AT-THE-END" BLUEPRINT
        Explain how to change the materials NOW so the item is safe in EVERY scenario:
        1. **If it goes to Waste:** How we change the chemicals/binders to ensure it safely mineralizes (Biological Nutrient).
        2. **If it goes to Recycle:** How we simplify the DNA (Mono-materials) to ensure it stays high-value.
        *Focus on removing 'Material Incompatibility' at the design stage.* (300+ Words).

        ### 🎯 THE BEST PATH & HANDLING
        Clearly identify the #1 most profitable path. Provide step-by-step handling and backhauling logistics for a company like Costco.

        ---

        ### 🏁 EXECUTIVE SUMMARY: THE SUSTAINABILITY LEADERSHIP VERDICT
        Explain why this design makes the company an industry leader. How does creating a 'consumer-proof' sustainable product eliminate future liability and maximize brand trust? (200 Words).
        """

        st.markdown(generate_vora_analysis(master_prompt))

# --- 5. FOOTER ---
st.sidebar.info("VoraCycle Enterprise Logic: v4.6.0")