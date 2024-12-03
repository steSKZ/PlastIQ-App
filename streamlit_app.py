import streamlit as st

# Set page configuration
st.set_page_config(page_title="PlastIQ - Result Page", layout="wide")

page_input_contact = st.Page("subpages/input/contact.py", title="Kontakt", icon=":material/contacts:")
page_input_company = st.Page("subpages/input/company.py", title="Unternehmen", icon=":material/apartment:")
page_input_product = st.Page("subpages/input/product.py", title="Zusammensetzung", icon=":material/add_circle:")
page_input_pr_origin = st.Page("subpages/input/product_origin.py", title="Herkunft", icon=":material/add_circle:")
page_input_pr_quality = st.Page("subpages/input/product_quality.py", title="Qualität – Wertstoff", icon=":material/add_circle:")
page_input_pr_quality_additive = st.Page("subpages/input/product_quality_additive.py", title="Qualität – Zusätze", icon=":material/add_circle:")
page_input_pr_further = st.Page("subpages/input/product_further.py", title="Weitere Angaben", icon=":material/add_circle:")
page_output_score = st.Page("subpages/output/score.py", title="Passende Verwertung", icon=":material/all_match:")


pg = st.navigation(
    {
        "Datenaufnahme - Kontakt": [page_input_contact, page_input_company],
        "Datenaufnahme - Wertstoff": [page_input_product, page_input_pr_origin, page_input_pr_quality, page_input_pr_quality_additive, page_input_pr_further],
        "Ergebnisse": [page_output_score],
    }
)

# Initialize keys for product input form if not available 
if "key_dict_product" not in st.session_state:
    st.session_state.key_dict_product = {"input_waste_fraction_number":1}

    for i in range(4):
        st.session_state.key_dict_product[f"input_wertstoff_typ_{i}"] = None
        st.session_state.key_dict_product[f"input_wertstoff_name_{i}"] = ""
        st.session_state.key_dict_product[f"input_wertstoff_anteil_{i}"] = 0.00
    

# Initialize keys for product input form if not available 
if "key_dict_product_origin" not in st.session_state:
    st.session_state.key_dict_product_origin = {"input_wertstoff_origin":None,
                                                "input_wertstoff_use":"",
                                                "input_wertstoff_collection":None,
                                                "input_wertstoff_code":"",
                                                }
    
# Initialize keys for product input form if not available 
if "key_dict_product_quality" not in st.session_state:
    st.session_state.key_dict_product_quality = {"input_wertstoff_reach":None,
                                                "input_wertstoff_colour":None,
                                                "input_wertstoff_purity":100,
                                                "input_wertstoff_contaminants_level":None,
                                                "input_wertstoff_contaminants_type":"",
                                                }

    for i in range(4):
        st.session_state.key_dict_product_quality[f"input_additiv_typ_{i}"] = []
        st.session_state.key_dict_product_quality[f"input_fuellstoff_typ_{i}"] = []

pg.run()