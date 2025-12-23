import streamlit as st
from openai import OpenAI

# --- 1. ACCESS & CONNECTION ---
# Fixed the KeyError by using a single-layer secret call
try:
    client = OpenAI(api_key=st.secrets["OPENAI_API_KEY"])
except Exception as e:
    st.error("🚨 API Key Missing: Please add 'OPENAI_API_KEY' to your Streamlit Secrets.")
    st.stop()

# --- 2. THE ARBITER ENGINE ---
def generate_vora_analysis(prompt):
    """Connects to GPT-4o with high-creativity temperature for strategic variety."""
    try:
        response = client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are the VoraCycle Lead Arbiter. You provide high-density, forensic-level strategic audits. You never give short, generic answers."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7 # Strategic variety sweet spot
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Analysis Error: {str(e)}"

# --- 3. USER INTERFACE ---
st.set_page_config(page_title="VoraCycle Arbiter", layout="wide")
st.title("🛡️ VoraCycle: Pre-Lifecycle Forensic Arbiter")
st.markdown("### Determining the Path of Maximum Capital & Environmental Retention")

query = st.text_input("Enter the Asset, Material, or Product for Strategic Audit:", 
                     placeholder="e.g., Industrial Lithium-Ion Battery Housing")

if query:
    st.divider()
    
    with st.spinner("⚖️ VoraCycle Arbiter is calculating the optimal circular path..."):
        # --- 4. THE HIGH-DENSITY PROMPT ---
        master_prompt = f"""
        Act as the VoraCycle Lead Arbiter. Conduct a Pre-Lifecycle Strategic Audit for: {query}.
        
        GOAL: Determine the strategy that maximizes BOTH Business Profit (ROI) and Environmental Restoration.
        
        REQUIRED OUTPUT FORMAT:
        
        ### 🎯 STRATEGIC SCORECARD
        - **Strategic Confidence Score:** [Provide 0-100%]
        - **Business Value Rating:** [High/Medium/Low]
        - **Environmental Recovery Potential:** [High/Medium/Low]
        - **Recommended Path:** [Refurbish, Repurpose, or (as a last resort) Mineralize]

        ---

        ### 🧪 THE 'WHY': ENVIRONMENTAL & BUSINESS SYNERGY
        Provide a 250-word deep-dive forensic analysis. Explain why this specific path is the 'Highest Form of Sustainability.' 
        Discuss 'Asset Lifecycle Extension,' 'Residual Value Recovery,' and 'Capital Resilience.' 
        Explain exactly how this strategy prevents 'Capital Leakage' while maintaining high operational velocity. 
        Focus on how the business wins by NOT turning this into waste.

        ### ⚖️ BEST STRATEGIC PATH: COMPETITIVE JUSTIFICATION
        Provide a 250-word technical comparison. Explain why you REJECTED other options (especially Mineralization). 
        Why is this path superior for the planet and the P&L? Discuss 'Carbon Substitution' and 
        'Supply Chain Sovereignty.' Prove that this path is the most mathematically sound decision 
        for a company seeking long-term stability.

        ### 🛡️ ENVIRONMENTAL FAILSAFE & 5-YEAR OUTCOME
        Provide a 150-word summary. If this asset eventually reaches end-of-life, how does it return 
        to the biological cycle safely? What is the projected 5-year systemic impact of choosing this path now?
        """

        # --- 5. EXECUTION & DISPLAY ---
        analysis_result = generate_vora_analysis(master_prompt)
        
        # Display the result in a clean, professional container
        st.markdown(analysis_result)

# --- 6. FOOTER ---
st.sidebar.info("VoraCycle Arbiter Logic: v2.5.0 (High-Density Forensic Mode)")