import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import pandas as pd
from PIL import Image
import io

# Page Configuration
st.set_page_config(
    page_title="Erythromycin Analysis",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS for Premium Look & Header Positioning
st.markdown("""
<style>
    .main {
        background-color: #0e1117;
        color: #ffffff;
    }
    .student-info {
        position: absolute;
        top: -60px;
        right: 10px;
        text-align: right;
        line-height: 1.2;
        font-size: 0.9em;
        opacity: 0.8;
    }
    .metric-card {
        background-color: #1e2130;
        padding: 20px;
        border-radius: 10px;
        border-left: 5px solid #ff4b4b;
        margin-bottom: 20px;
    }
    .centered-title {
        text-align: center;
        margin-top: 20px;
        margin-bottom: 40px;
    }
</style>

<div class="student-info">
    <b>Name:</b> [YOUR NAME]<br>
    <b>Class/Section:</b> [CLASS/SECTION]<br>
    <b>Roll No:</b> [ROLL NO.]
</div>
""", unsafe_allow_html=True)

# Main Title (Top Middle)
st.markdown("<h1 class='centered-title'>Erythromycin Chiral Center Analysis</h1>", unsafe_allow_html=True)

# Hardcoded Erythromycin SMILES
ERYTHROMYCIN_SMILES = "CC[C@H]1[C@@H](C)[C@@H](O)[C@H](OC(=O)[C@H](C)[C@@H](O)[C@H](C)O)[C@@H](C)[C@@H](O)[C@H](C)O[C@H]2[C@H](O)[C@@H](C)O[C@@H](O[C@@H]3[C@H](C)[C@@H](O)[C@H](C)O[C@H]3O)[C@H](O)[C@H]2O[C@@H]1C"

try:
    # Create molecule
    mol = Chem.MolFromSmiles(ERYTHROMYCIN_SMILES)
    
    if mol:
        # 1. Molecule Visualization & Summary
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.subheader("Chemical Structure")
            img = Draw.MolToImage(mol, size=(600, 400), kekulize=True)
            st.image(img, use_container_width=True, caption="2D Structure of Erythromycin")

        # 2. Chiral Center Calculation
        mol_with_hs = Chem.AddHs(mol)
        Chem.AssignStereochemistry(mol_with_hs, force=True, cleanIt=True)
        chiral_centers = Chem.FindMolChiralCenters(
            mol_with_hs,
            includeUnassigned=False,
            useLegacyImplementation=False
        )
        
        with col2:
            st.subheader("Statistical Summary")
            st.markdown(f"""
            <div class="metric-card">
                <h2 style='margin:0; color:#ff4b4b;'>{len(chiral_centers)}</h2>
                <p style='margin:0; opacity:0.8;'>Total Chiral Centers Detected</p>
            </div>
            """, unsafe_allow_html=True)
            
            # Basic Properties
            mw = Chem.rdMolDescriptors.CalcExactMolWt(mol)
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            
            st.write(f"**Formula:** {formula}")
            st.write(f"**Molecular Weight:** {mw:.2f} g/mol")

        st.markdown("---")
        
        # 3. Detailed Results
        st.subheader("Detailed Chiral Center Mapping")
        if len(chiral_centers) > 0:
            data = []
            for idx, config in chiral_centers:
                atom = mol_with_hs.GetAtomWithIdx(idx)
                element = atom.GetSymbol()
                data.append({
                    "Atom Index": idx,
                    "Atom Element": element,
                    "Stereo Configuration (R/S)": config
                })
            
            df = pd.DataFrame(data)
            st.table(df)

except Exception as e:
    st.error(f"Error processing molecule details: {str(e)}")

# Footer
st.markdown("---")
st.caption("Chemistry Project | Powered by RDKit & Streamlit")

# Footer
st.markdown("---")
st.caption("Powered by RDKit & Streamlit | Chemistry Project v2")
