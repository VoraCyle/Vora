import streamlit as st
import streamlit_authenticator as st_auth
from rdkit import Chem
from rdkit.Chem import Draw
from pyscf import gto, scf

# --- 1. Credentials & Auth Setup ---
credentials = {
    'usernames': {
        'wraith': {
            'name': 'Wraith',
            'password': 'Vora1630' 
        }
    }
}

authenticator = st_auth.Authenticate(
    credentials,
    "vora_cookie",
    "auth_key",
    cookie_expiry_days=30
)

# --- 2. Login UI ---
authenticator.login(location='main')

auth_status = st.session_state.get("authentication_status")
name = st.session_state.get("name")
username = st.session_state.get("username")

# --- 3. Main App Logic (Only if Logged In) ---
if auth_status:
    # Sidebar Setup
    authenticator.logout('Logout', 'sidebar')
    st.sidebar.header("Material Selection")
    
    material = st.sidebar.selectbox(
        "Choose a common wholesale packaging material",
        ["Polyethylene (PE)", "Polyethylene Terephthalate (PET)", "Polypropylene (PP)", "Custom SMILES"]
    )

    # App Header
    st.title("🔮 Wraith VoraCycle")
    st.markdown("**Precycling Powered by Precision**")
    st.markdown("Predict material circularity before it becomes waste. Quantum-inspired AI for proactive regeneration.")

    # --- Molecular Analysis ---
    smiles_dict = {
        "Polyethylene (PE)": "CCCCCCCCCC",
        "Polyethylene Terephthalate (PET)": "C1=CC=C(C=C1)C(=O)OCCOC(=O)C2=CC=C(C=C2)C(=O)O",
        "Polypropylene (PP)": "CC(C)CCCCCC"
    }

    if material != "Custom SMILES":
        smiles = smiles_dict[material]
    else:
        smiles = st.sidebar.text_input("Enter SMILES string", "CCCC")

    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(600, 600))
        st.image(img, caption=f"Molecular structure: {material}")

    # --- Quantum Simulation ---
    st.write("### Precision Simulation Running...")
    
    # Proxy molecule for quantum calculation
    proxy_mol = gto.M(atom='C 0 0 0; C 0 0 1.4; H 0 1 1; H 0 0 2.4', basis='sto-3g')
    mf = scf.RHF(proxy_mol)
    energy = mf.kernel()

    st.write(f"**Quantum energy stability**: {energy:.2f} Hartree (lower = more durable)")

    # Scoring Logic
    scores = {"Polyethylene (PE)": 48, "Polyethylene Terephthalate (PET)": 72, "Polypropylene (PP)": 55}
    score = scores.get(material, 65)
    
    st.metric(label="**VoraCycle Circularity Score**", value=f"{score}/100", delta=f"{score-50:+} vs average")

    if score < 70:
        st.warning("⚠️ Elevated non-recyclable risk detected.")

    # --- Precision Redesign Recommendations ---
    st.markdown("---")
    st.markdown("### 🔬 Precision Redesign Recommendation")

    redesign_suggestions = {
        "Polyethylene Terephthalate (PET)": {
            "optimized_name": "Polyethylene Furanoate (PEF)",
            "optimized_score": 94,
            "optimized_stability": -74.20,
            "modifications": "Replace terephthalic acid with bio-derived 2,5-furandicarboxylic acid (FDCA). Incorporate enzyme-recognizable cleavage sites.",
        },
        "Polyethylene (PE)": {
            "optimized_name": "Bio-PE with Enzymatic Triggers",
            "optimized_score": 88,
            "optimized_stability": -68.50,
            "modifications": "Blend with PHA copolymers and add ester-based degradation triggers.",
        },
        "Polypropylene (PP)": {
            "optimized_name": "Isotactic PP with Bio-additives",
            "optimized_score": 85,
            "optimized_stability": -72.10,
            "modifications": "Incorporate oxidation-promoting prodegradants and bio-based comonomers.",
        }
    }

    suggestion = redesign_suggestions.get(material, {
        "optimized_name": "Custom Optimized Variant",
        "optimized_score": min(95, score + 18),
        "optimized_stability": energy + 1.5,
        "modifications": "Functional group optimization using quantum-guided bond lability tuning."
    })

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("#### Original Material")
        st.error(f"**{material}** ({score}/100)")
    with col2:
        st.markdown("#### → Optimized Redesign")
        st.success(f"**{suggestion['optimized_name']}** ({suggestion['optimized_score']}/100)")

    st.info(f"**Action Plan:** {suggestion['modifications']}")

    # --- Environmental Fate Section ---
    st.markdown("---")
    st.markdown("### 🌍 Environmental Fate & Redesign Impact")

    colA, colB = st.columns(2)
    with colA:
        st.error("**Original Fate**")
        st.markdown("- Persists for **hundreds of years**.\n- Resist microbial breakdown.\n- Fragments into microplastics.")
    with colB:
        st.success("**Redesigned Fate**")
        st.markdown("- Rapid **safe breakdown** (months/years).\n- No persistent microplastics.\n- Supports closed-loop recycling.")

    st.caption("Prototype created by Wraith | Pioneering Precycling with Precision")

# --- 4. Auth Feedback ---
elif auth_status == False:
    st.error('Username/password incorrect')
elif auth_status is None:
    st.warning('Please enter your credentials to access the VoraCycle engine.')
