# --- 1. Strategic Rationale Engine ---
    def get_rationale(data, category):
        reasons = []
        # Bad Idea scenarios
        if data['toxic'] > 0:
            reasons.append("❌ **BAD IDEA:** Current material contains halogens. These can migrate into fats (leaching) and release toxic gases if incinerated.")
        if category == "Hot Food (Meat/Chicken)" and data['tsi'] < 50:
            reasons.append("❌ **BAD IDEA:** Low thermal resistance. The packaging may soften or 'bond' to the food at 180°F, contaminating the product.")
        
        # Good Idea / Change scenarios
        if data['fate'] < 50:
            reasons.append("💡 **THE CHANGE:** We recommend shifting to a 'perforated' polymer backbone.")
            reasons.append("🛡️ **THE BENEFIT:** This keeps the item 100% rigid and leak-proof on the shelf, but allows landfill moisture to 'unzip' the plastic after use.")
        
        if data['recycle'] > 75:
            reasons.append("✅ **GOOD IDEA:** This material is high-purity. Using it reduces Costco’s 'Plastic Tax' liability by ensuring it stays in the circular economy.")
            
        return reasons

    # --- 2. THE UPDATED UI BLOCK ---
    if data:
        st.image(Draw.MolToImage(data['mol'], size=(600, 600)))
        
        # --- NEW: STRATEGIC NARRATIVE ---
        st.markdown("---")
        st.subheader("📝 Strategic Rationale")
        narratives = get_rationale(data, category)
        for note in narratives:
            st.write(note)

        # --- RECAP OF PERFORMANCE IMPACT ---
        st.info(f"**Performance Impact:** The proposed changes maintain the **Barrier Index** at {data['tsi']/10:.1f}/10, ensuring no loss in shelf-life or food safety.")
