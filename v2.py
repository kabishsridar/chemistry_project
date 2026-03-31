import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import pandas as pd
from PIL import Image
import io

# Page Configuration
st.set_page_config(
    page_title="Chiral center Analyzer",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for Premium Look
st.markdown("""
<style>
    .main {
        background-color: #0e1117;
        color: #ffffff;
    }
    .stButton>button {
        width: 100%;
        border-radius: 5px;
        height: 3em;
        background-color: #ff4b4b;
        color: white;
        font-weight: bold;
        border: none;
        transition: 0.3s;
    }
    .stButton>button:hover {
        background-color: #ff6b6b;
        box-shadow: 0 4px 15px rgba(255, 75, 75, 0.4);
    }
    .metric-card {
        background-color: #1e2130;
        padding: 20px;
        border-radius: 10px;
        border-left: 5px solid #ff4b4b;
        margin-bottom: 20px;
    }
    .molecule-container {
        background-color: white;
        padding: 20px;
        border-radius: 15px;
        display: flex;
        justify-content: center;
        margin-top: 20px;
    }
</style>
""", unsafe_allow_html=True)

# Sidebar
with st.sidebar:
    st.title("🧪 ChemAnalyzer")
    st.markdown("---")
    st.markdown("### Input Molecule")
    
    # Common molecules for quick selection
    presets = {
        "Erythromycin": "CC[C@H]1[C@@H](C)[C@@H](O)[C@H](OC(=O)[C@H](C)[C@@H](O)[C@H](C)O)[C@@H](C)[C@@H](O)[C@H](C)O[C@H]2[C@H](O)[C@@H](C)O[C@@H](O[C@@H]3[C@H](C)[C@@H](O)[C@H](C)O[C@H]3O)[C@H](O)[C@H]2O[C@@H]1C",
        "Lactic Acid": "CC(O)C(=O)O",
        "Glucose": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
    }
    
    selection = st.selectbox("Choose a preset:", list(presets.keys()))
    smiles_input = st.text_area("Or Enter SMILES String:", presets[selection], height=100)
    
    analyze_btn = st.button("Analyze Molecule")

# Main Content
st.title("Molecular Stereochemistry Analysis")
st.markdown("Analyze chiral centers and molecular properties instantly.")

if smiles_input:
    try:
        # Create molecule
        mol = Chem.MolFromSmiles(smiles_input)
        
        if mol is None:
            st.error("Invalid SMILES string. Please check your input.")
        else:
            # 1. Molecule Visualization
            col1, col2 = st.columns([1, 1])
            
            with col1:
                st.subheader("Structure")
                # Generate image
                img = Draw.MolToImage(mol, size=(600, 400), kekulize=True)
                st.image(img, use_container_width=True, caption=f"2D Representation of {selection if smiles_input == presets[selection] else 'Input'}")

            # 2. Chiral Center Calculation
            mol_with_hs = Chem.AddHs(mol)
            Chem.AssignStereochemistry(mol_with_hs, force=True, cleanIt=True)
            chiral_centers = Chem.FindMolChiralCenters(
                mol_with_hs,
                includeUnassigned=False,
                useLegacyImplementation=False
            )
            
            with col2:
                st.subheader("Summary")
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
            st.subheader("Chiral Centers Details")
            if len(chiral_centers) > 0:
                data = []
                for idx, config in chiral_centers:
                    atom = mol_with_hs.GetAtomWithIdx(idx)
                    element = atom.GetSymbol()
                    data.append({
                        "Atom Index": idx,
                        "Element": element,
                        "Configuration (R/S)": config
                    })
                
                df = pd.DataFrame(data)
                st.table(df)
            else:
                st.info("No chiral centers found in this molecule.")

    except Exception as e:
        st.error(f"Error processing molecule: {str(e)}")
else:
    st.info("Please enter a SMILES string in the sidebar to begin.")

# Footer
st.markdown("---")
st.caption("Powered by RDKit & Streamlit | Chemistry Project v2")
