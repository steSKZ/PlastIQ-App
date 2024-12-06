## Import libraries
import streamlit as st
import pandas as pd

## Define variables
# Initialize variables
wertstoffscore = 0 # variable for the final score of the waste material
percent_valuable_substance = 0.0 #set up variable for the amount of valuable substance in the waste
material_type = [] # list of waste types added during the input
material_share = [] # list of corresponding waste percentages added during the input

# Parameters changing model
set_threshold_assessibility = 0.8 # threshold for Step 2 regarding the assessibility of the recycling fraction
wt_ferromagnetic = 1 # weigth of sorting method "ferromagnetic", between 0 and 1
wt_eddycurrent = 1 # weigth of sorting method "eddycurrent", between 0 and 1
wt_density = 1 # weigth of sorting method "density", between 0 and 1
wt_electrostatic = 1 # weigth of sorting method "electrostatic", between 0 and 1
wt_sensorbased = 1 # weigth of sorting method "sensorbased", between 0 and 1

# Default LCA data TODO change to tuple
lca_declared_unit = 1 # Declared unit, usually 1 kg
lca_emission_origen = "DE"

lca_substitution_factor_plastic = 0.5 # substitution factor for avoided emissions for plastics # TODO 
lca_substitution_factor_metal = 0.8 # substitution factor for avoided emissions for metal # TODO 
lca_substitution_factor_electricity = 1 # substitution factor for avoided emissions for electricity # TODO
lca_substitution_factor_heat = 1 # substitution factor for avoided emissions for heat # TODO 

lca_distance_to_wte = 10 # Distance to waste-to-energy-plant in km TODO
lca_distance_to_landfill = 10 # Distance from wte-plant to landfill in km TODO
lca_distance_to_recycling_metal = 10 # Distance from wte-plant to recycling facility for metal in km TODO

lca_vehicle_to_wte = "transport_lkw_22" # Chosen vehicle to transport materials to wte
lca_vehicle_to_landfill = "transport_lkw_22" # Chosen vehicle to transport materials to landfill
lca_vehicle_to_recycler_metal = "transport_lkw_22" # Chosen vehicle to transport materials to recycler for metal
lca_vehicle_to_recycler_plastic = "transport_lkw_22" # Chosen vehicle to transport materials to recycler for plastic

lca_efficiency_wte_electric = 0.113
lca_efficiency_wte_heat = 0.33
lca_share_dross = 0.2 # TODO share of dross in ratio to 1 kg of burnable material (e.g. plastic)

lca_electricity_use_sorting_metal = 0.3 # electricity use for sorting metal per kg, in kWh 
lca_electricity_use_shredding_metal = 0.3 # electricity use for shredding metal per kg, in kWh TODO
lca_electricity_use_cleaning_metal = 0.3 # electricity use for cleaning metal per kg, in kWh TODO
lca_electricity_use_melting_metal = 1 # electricity use for melting metal per kg, in kWh TODO
lca_heat_use_cleaning_metal = 1 # heat energy use for cleaning metal per kg, in kWh TODO
lca_water_use_cleaning_metal = 1 # water use for cleaning metal per kg in kg Water TODO
lca_wastewater_use_cleaning_metal = 1 # wastewater for cleaning metal per kg in kg wastewater TODO
lca_solvent_use_cleaning_metal = 1 # solvent use for clening metal per kg in kg solvent TODO

lca_electricity_use_sorting_mixed = 0.3 # electricity use for sorting mixed materials per kg in kWh TODO
lca_electricity_use_shredding_plastic = 0.3 # electricity use for shredding plastic per kg in kWh TODO
lca_electricity_use_cleaning_plastic = 0.3 # electricity use for cleaning plastic per kg in kWh TODO
lca_heat_use_cleaning_plastic = 1 # heat energy use for cleaning metal per kg, in kWh TODO
lca_water_use_cleaning_plastic = 1 # water use for cleaning metal per kg in kg Water TODO
lca_wastewater_use_cleaning_plastic = 1 # wastewater for cleaning metal per kg in kg wastewater TODO
lca_solvent_use_cleaning_plastic = 1 # solvent use for clening metal per kg in kg solvent TODO

lca_share_io_shredding = 0.99 # share of i.O. material fit for recycling after the shredding process TODO
lca_share_io_cleaning = 0.99 # share of i.O. material fit for recycling after the cleaning process TODO
lca_share_io_regranulation = 0.99 # share of i.O. material fit for recycling after the regranulation process TODO
lca_share_io_overall = lca_share_io_shredding * lca_share_io_cleaning * lca_share_io_regranulation # share of i.O. material fit for recycling overall

# file paths
file_path_background = "content/background_data_decision_tree.xlsx"
file_path_input = "content/plastiq_input_information.xlsx"

# Define labels for different dataframes
# for dataframe materials
label_material_type = "type"
label_material_category = "category"
label_material_share = "share"
label_material_result_sorting = "result_sorting"
columns_materials = [label_material_type, label_material_category, label_material_share, label_material_result_sorting]

# for dataframe of sorting results:
label_material_pairing = "material_pairing"
label_sort_ferromagnetic = "sort_ferromagnetic"
label_sort_eddycurrent = "sort_eddycurrent"
label_sort_density = "sort_density"
label_sort_electrostatic = "sort_electrostatic"
label_sort_sensorbased = "sort_sensorbased"
columns_sorting = [label_material_pairing, label_sort_ferromagnetic, label_sort_eddycurrent, label_sort_density, label_sort_electrostatic, label_sort_sensorbased]

# for dataframe emissions
columns_emissions = [
    "emissions_overall_per_year", "emissions_overall_per_weight", "emissions_overall_per_declared_unit",
    "processing_plastic_total", "processing_plastic_sorting", "processing_plastic_shredding", "processing_plastic_cleaning", "processing_plastic_drying", "processing_plastic_regranulation", 
    "processing_metal_total", "processing_metal_sorting", "processing_metal_shredding", "processing_metal_cleaning", "processing_metal_melting",
    "endoflife_total", "endoflife_incineration", "endoflife_landfill",
    "transport_total", "transport_recycler_plastic", "transport_recycler_metal", "transport_wte", "transport_landfill",
    "avoided_total", "avoided_electricity", "avoided_heat", "avoided_plastic", "avoided_metal"
] 

## Get session state variables
# Get weight of materials
waste_weigth = st.session_state.key_dict_product_amount["input_wertstoff_menge"]
# Get weight of material per year
material_weight_regular = st.session_state.key_dict_product_amount["input_haeufigkeit_menge"]
material_frequency = st.session_state.key_dict_product_amount["input_haeufigkeit_turnus"]
# Get reach conformity
reach_conformity = st.session_state.key_dict_product_quality["input_wertstoff_reach"]
# Get number of fractions
material_fraction_number = st.session_state.key_dict_product["input_waste_fraction_number"]
# Get material type and share of fraction to total material
for k in range(material_fraction_number):
    material_type.append(st.session_state.key_dict_product[f"input_wertstoff_typ_{k}"])
    material_share.append(st.session_state.key_dict_product[f"input_wertstoff_anteil_{k}"])

## Get input data from files
# Get list with all available materials and data on them
df_materials = pd.read_excel(file_path_background, sheet_name = "list_material")
# Magnetabscheidung
df_sort_ferromagnetic = pd.read_excel(file_path_background, sheet_name = "sort_ferromagnetic", index_col=0)
# Wirbelstromsortierung
df_sort_eddycurrent = pd.read_excel(file_path_background, sheet_name = "sort_eddycurrent", index_col=0)
# Dichtesortierung
df_sort_density = pd.read_excel(file_path_background, sheet_name = "sort_density", index_col=0)
# Elektrostatische Sortierung
df_sort_electrostatic = pd.read_excel(file_path_background, sheet_name = "sort_electrostatic", index_col=0)
# Sensorbased Sorting 
df_sort_sensorbased = pd.read_excel(file_path_background, sheet_name = "sort_sensorbased", index_col=0)
# Get emission data
df_lca_emission_data = pd.read_excel(file_path_background, sheet_name = "lca_calculation")

## Define functions - general
# Function to locate WS to a class
def func_evaluateWS(wertstoffscore: float):
    ws_rounded = round(wertstoffscore, 1)
    step = 0.1
    if -10 <= ws_rounded < 0: #Klasse G: Sondermüll
      ws_category = "H"
      ws_sentence = "Dein Abfall sollte aufgrund der Inhalte auf dem Sondermüll entsorgt werden."
    elif 0 <= ws_rounded < 50: #Klasse F: Energetische Verwertung
      ws_category = "F"
      ws_sentence = "Dein Abfall eigent sich voraussichtlich nur für die energetische Verwertung."
    elif 50 <= ws_rounded < 60: #Klasse E: Recycling, chemisch
      ws_category = "E"
      ws_sentence = "Deine Wertstoffe eigenen sich voraussichtlich für das chemische Recycling."
    elif 60 <= ws_rounded < 75: #Klasse D: Recycling, werkstofflich, gemischt
      ws_category = "D"
      ws_sentence = "Deine Wertstoffe eignen sich voraussichtlich für das gemischte mechanische Recycling."
    elif 75 <= ws_rounded < 90: #Klasse C: Recyling, werkstofflich, sortenrein
      ws_category = "C"
      ws_sentence = "Deine Wertstoffe eignen sich voraussichtlich für das sortenreine mechanische Recycling."
    elif 90 <= ws_rounded < 95: #Klasse B: Wiederverwendung, niederwertig
      ws_category = "B"
      ws_sentence = "Deine Wertstoffe können mit einer anderen Funktion wiederverwertet werden."
    elif 95 <= ws_rounded <= 100: #Klasse A: Wiederverwendung, gleichwertig
      ws_category = "A"
      ws_sentence = "Deine Wertstoffe können mit der gleichen Funktion wiederverwertet werden."
    else: #Fehlerhafte Eingabe
      ws_category = "Eingabe konnte nicht verarbeitet werden."
    return ws_category, ws_sentence

# Initialize dataframe for sorting output scores
def func_initializeOutputDf(amount: float = material_fraction_number, label_columns: list = columns_sorting) -> pd.DataFrame:
    df_new = pd.DataFrame(columns=label_columns)
    list_material_pairing = []

    # Loop through each entry and compare it with every other entry, while avoiding duplicates
    # Add row for each material pairing
    for k in range(amount):
        for l in range(k+1, amount):
            
            # Specify the row and column indices
            first_material = material_type[k]
            second_material = material_type[l]

            # Define material pairing
            pair_string = f"{first_material}_{second_material}"

            # Add to list
            list_material_pairing.append(pair_string)

    # Add list to dataframe
    df_new[label_material_pairing] = list_material_pairing

    # Output new dataframe
    return df_new

# Initialize dataframe for result 
def func_initializeResultDf(amount: float = material_fraction_number, label_columns: list = columns_materials) -> pd.DataFrame:
    
    # New dataframe with material as column name
    df_new = pd.DataFrame(columns = label_columns)
    # Add first column "material type"
    df_new[label_material_type] = material_type
    # Lookup category of material
    material_category = []

    for k in range(amount):
        
        # lookup category from material dataframe
        category = df_materials.loc[df_materials["abbreviation"].str.fullmatch(material_type[k], case=False, na=False), "category"]
        material_category.append(category.iloc[0])

    # Add second column with material category
    df_new[label_material_category] = material_category

    # Add third column with material share
    df_new[label_material_share] = material_share            

    return df_new

def func_get_indices_for_material(df: pd.DataFrame, column: str, material: str) -> list:
    """
    Get the indices of rows where a specific material is present in a specified column.

    Parameters:
        df (pd.DataFrame): The DataFrame to search.
        column (str): The column name to search in.
        material (str): The material to look for.

    Returns:
        list: A list of row indices where the material is present.
    """
    indices = df[df[column].str.contains(material, case=False, na=False)].index
    return indices.tolist()

# Function to check all available materials for their sorting options
def func_checkSorting(df: pd.DataFrame, amount: float = material_fraction_number) -> list:
    
    # Initialize list with all result values
    values_sort = []

    # Loop through each entry and compare it with every other entry, while avoiding duplicates
    for k in range(amount):
        for l in range(k+1, amount):
              
          # Specify the row and column indices
          row_label = material_type[k]
          column_label = material_type[l]

          # Get and return value at indices from matrix
          result = df.loc[row_label, column_label]

          # Add to list
          values_sort.append(float(result))
    
    # Output list with result values
    return values_sort

def func_sum_columns_with_conditions(df, row_label, include_substring, exclude_substring) -> float:
    """
    Sums values in a row for columns containing a specific substring
    but excluding columns containing another substring.

    Parameters:
        df (pd.DataFrame): The DataFrame to process.
        row_label (str): The label of the row to sum.
        include_substring (str): Substring to include in column names.
        exclude_substring (str): Substring to exclude from column names.

    Returns:
        float: The sum of the selected values.
    """
    filtered_columns = [
        col for col in df.columns
        if include_substring in col and exclude_substring not in col
    ]
    return df.loc[row_label, filtered_columns].sum()

## Define functions - LCA
# Get a specific emission factor from background_data
def func_lca_get_emission_factor(label_category: str, label_use: str) -> float:
    """
    Get a specific emission factor from the background data based on category and use label.

    Parameters:
        label_category (str): The category to filter the emission data by.
        label_use (str): The use label to search for in the filtered emission data.

    Returns:
        float: The emission factor (GWP100) corresponding to the specified category and use label.

    Raises:
        ValueError: If no emission factor is found for the given category and use label.
    """
    # Filter the DataFrame by the specified category
    df_filtered = df_lca_emission_data[df_lca_emission_data["category"] == label_category]

    # Lookup emission factor specific to the label_use
    emission_factor = df_filtered.loc[
        df_filtered["use for"].str.contains(label_use, case=False, na=False), "GWP100"
    ]

    # Check if an emission factor was found
    if emission_factor.empty:
        raise ValueError(
            f"No emission factor found for category '{label_category}' and use label '{label_use}'."
        )

    # Return the first matching emission factor (or handle multiple matches if required)
    return float(emission_factor.iloc[0])

# get combined share of a specific category (plastic, meta)
def func_lca_get_category_share(df: pd.DataFrame, relevant_category: str) -> tuple:
    """
    Calculates the combined share of a specific category in a DataFrame and returns
    both the filtered DataFrame and the combined share.

    Parameters:
        df (pd.DataFrame): The DataFrame containing material data.
        relevant_category (str): The category for which the share is calculated (e.g., 'plastic', 'metal').

    Returns:
        tuple: A tuple containing:
            - float: The total combined share of the specified category.
            - pd.DataFrame: The filtered DataFrame containing only rows of the relevant category.

    Raises:
        ValueError: If the DataFrame does not contain the required columns or the category is not found.
    """
    # Ensure required columns exist in the DataFrame
    required_columns = [label_material_category, label_material_share]
    if not all(col in df.columns for col in required_columns):
        raise ValueError(f"The DataFrame must contain the columns: {required_columns}")

    # Filter DataFrame for the relevant category
    df_filtered = df[df[label_material_category] == relevant_category]

    # Sum of the share of relevant materials
    category_share = df_filtered[label_material_share].sum()

    # Return the combined share and the filtered DataFrame
    return float(category_share/100), df_filtered

# Get emissions from transport
def func_lca_emissions_transport(label_use: str, distance: float = 100.0, payload: float = 1.0) -> float:
    # get emission factor
    emission_factor = func_lca_get_emission_factor(label_category="transport", label_use=label_use)
    # calculate emissions with payload (in ton) * distance (in km) * emission factor (in kg CO2e/(t*km))
    emission_transport = payload * distance * emission_factor
    
    return emission_transport

# get emissions from waste incineration
def func_lca_emissions_incineration(df: pd.DataFrame, weight: float) -> tuple:
    
    # Initialize values for resulting emissions
    emissions_incineration = 0 #in kg CO2e/kg waste
    
    # for every material fraction
    for k in range(len(df[label_material_type].tolist())):
        
        # get material name and percentage from df_result in a list [name, percentage]
        material_type = df.at[k, label_material_type]
        material_share = float(df.at[k, label_material_share])

        # lookup name in df_emission, get emission_factor and add to list
        emission_factor = func_lca_get_emission_factor("incineration", material_type)

        # add to exisiting emissions and energy
        emissions_incineration += weight * material_share/100 * emission_factor
  
    return emissions_incineration

# get energy from waste incineration
def func_lca_energy_incineration(df: pd.DataFrame, weight: float) -> tuple:
    
    # Initialize values for electric and heat energy
    electric_energy_incineration = 0 #in kWh 
    heat_energy_incineration = 0 #in kWh
    
    # for every material fraction
    for k in range(len(df[label_material_type].tolist())):
        
        # get material name and percentage from df_result in a list [name, percentage]
        material_type = df.at[k, label_material_type]
        material_share = float(df.at[k, label_material_share])

        # lookup heating value from df_material
        lower_heating_value = float(df_materials.loc[df_materials["abbreviation"].str.fullmatch(material_type, case=False, na=False), "lower_heating_value_MJ_per_kg"])

        # calculate electric and heat energy [kWh] from heating value [J/kg], weight [kg], material_share [-] and efficiency [-]
        electric_energy_incineration += lower_heating_value * weight * material_share/100 * lca_efficiency_wte_electric / 3.6
        heat_energy_incineration += lower_heating_value * weight * material_share/100 * lca_efficiency_wte_heat / 3.6

    return electric_energy_incineration, heat_energy_incineration

# get emissions from process (categories: electricity, heat, ressources)
def func_lca_emissions_process(category: str, origen: str, source: str, amount: float = 1.0) -> float:
    """
    Calculates emissions from a process based on the category, origin, source, 
    and the amount of the process.

    Parameters:
        category (str): The category of the process (e.g., 'electricity', 'heat', 'resources').
        origen (str): The origin of the emission factor (e.g., geographical location).
        source (str): The specific source of emissions (e.g., type of fuel or material).
        amount (float): The process amount in relevant units (default is 1.0).

    Returns:
        float: The total emissions from the process.

    Raises:
        ValueError: If no emission factor is found for the specified category, origin, or source.
    """
    try:
        # Define the label for emission factor search
        label_use = f"{category}_{origen}_{source}"

        # Look up emission factor
        emission_factor = func_lca_get_emission_factor(label_category=category, label_use=label_use)

        # Calculate emissions from process amount and emission factor
        emission_process = amount * emission_factor

        return emission_process
    except ValueError as e:
        # Provide detailed error message if emission factor lookup fails
        raise ValueError(
            f"Failed to calculate emissions for category '{category}', origin '{origen}', "
            f"source '{source}', and amount {amount}. Error: {e}"
        )

# get cumulated emissions for plastic processes (drying, regranulation)
def func_lca_emissions_processes_plastic(process: str, df: pd.DataFrame, source_electricity: str = "mix", source_heat: str = "mix") -> float:

    # Initialize values for process emissions in kg CO2e per kWh
    emission_process_electricity = 0 
    emission_process_heat = 0 
    relevant_category = "plastic"

    # filter dataframe for relevant category
    category_share, df_filtered = func_lca_get_category_share(df=df, relevant_category=relevant_category)
        
    # set labels depending on process
    match process:
        case "drying":
            column_label_electricity = "drying_electricity_use_kWh_per_kg"
            column_label_heat = "drying_heat_use_kWh_per_kg"

        case "regranulation":
            column_label_electricity = "regranulation_electricity_use_kWh_per_kg"
            column_label_heat = "regranulation_heat_use_kWh_per_kg"

    # for every material fraction
    for k in range(df_filtered.shape[0]):
        
        # get material name and percentage from df_result in a list [name, percentage]
        material_type = df_filtered.iloc[k][label_material_type]
        material_share = df_filtered.iloc[k][label_material_share]

        # lookup electricity use for each process from list_material
        electricity_use = float(df_materials.loc[df_materials["abbreviation"].str.fullmatch(material_type, case=False, na=False), column_label_electricity])
        heat_use = float(df_materials.loc[df_materials["abbreviation"].str.fullmatch(material_type, case=False, na=False), column_label_heat])

        # calculate emissions [kg CO2e/kg/share plastic] for the relevant process from energy use [kWh], material share [-] and emission factor [kg CO2/kWh]
        emission_process_electricity += electricity_use * (material_share/category_share) * func_lca_emissions_process(category="electricity", origen=lca_emission_origen, source=source_electricity)
        emission_process_heat += heat_use * (material_share/category_share) * func_lca_emissions_process(category="heat", origen=lca_emission_origen, source=source_heat)

    # add up both emission sources
    emission_process = emission_process_electricity + emission_process_heat

    return emission_process

# get cumulated emissions for melting process
def func_lca_emissions_process_melting(df: pd.DataFrame) -> float:
    
    # initialize variable for process emissions
    emission_process = 0
    relevant_category = "metal"

    # filter dataframe for relevant category
    category_share, df_filtered = func_lca_get_category_share(df=df, relevant_category=relevant_category)
    
    # loop through every type of material
    for k in range(df_filtered.shape[0]):

        # get material type, category and share from df_result
        material_type = df_filtered.iloc[k][label_material_type]
        material_share = df_filtered.iloc[k][label_material_share]

        # calculate and add process emission of material
        emission_process +=  (material_share/category_share) * func_lca_emissions_process(category=f"melting-{relevant_category}", origen=lca_emission_origen, source=material_type)

    return emission_process

# get avoided emissions for secondary materials 1. plastic, 2. metal, depending on material in stream
def func_lca_emissions_avoided_material(df: pd.DataFrame, relevant_category: str) -> float:
    
    # initialize variable for avoided emissions
    emission_avoided_material = 0

    # filter dataframe for relevant category
    category_share, df_filtered = func_lca_get_category_share(df=df, relevant_category=relevant_category)

    # loop through every type of material
    for k in range(df_filtered.shape[0]):

        # get material type, category and share from df_result
        material_type = df_filtered.iloc[k][label_material_type]
        material_share = df_filtered.iloc[k][label_material_share]

        # calculate and add avoided emission of material
        emission_avoided_material +=  (material_share/category_share) * func_lca_emissions_process(category=f"production-{relevant_category}", origen=lca_emission_origen, source=material_type)

    return emission_avoided_material

# calculate material weight per year in ton
def func_lca_get_weigth_per_year(weight: float, frequency: str):
    # differentiate in calculation depending on the frequency statement
    match frequency:
      case "Tag":
          weight_per_year = weight * 365
      case "Woche":
          weight_per_year = weight * 52
      case "Monat":
          weight_per_year = weight * 12
      case "Quartal":
          weight_per_year = weight * 4
      case  "Jahr":
          weight_per_year = weight
      case _:
          weight_per_year = 0 #TODO Error Message
    return weight_per_year

# function to calculate emissions for scenario: materials go to waste-to-energy-plant
def func_lca_emissions_scenario_wte(df: pd.DataFrame, df_emission: pd.DataFrame) -> pd.DataFrame: 
    
    # determine the new index to add the scenario as new row to the emission dataframe
    new_index = len(df_emission)

    # transport to waste-to-energy(wte)-plant
    df_emission.at[new_index, "transport_wte"] = func_lca_emissions_transport(lca_vehicle_to_wte, lca_distance_to_wte, lca_declared_unit)

    # incineration of waste in waste-to-energy(wte)-plant
    df_emission.at[new_index, "endoflife_incineration"] = func_lca_emissions_incineration(df, lca_declared_unit)
    energy_incineration = func_lca_energy_incineration(df, lca_declared_unit)

    # sorting process after incineration between metals (recycling) and dross (landfill)
    share_metal, _ = func_lca_get_category_share(df_result, "metal") # share for remaining metal compared to 1 kg material
    share_plastic, _ = func_lca_get_category_share(df_result, "plastic") # share for remaining plastic compared to 1 kg material
    share_dross = lca_share_dross * share_plastic # share for the remaining dross from burned plastic compared to 1 kg waste

    # transport of dross to landfill + enmissions from landfill
    df_emission.at[new_index, "transport_landfill"] = share_dross * func_lca_emissions_transport(lca_vehicle_to_landfill, lca_distance_to_landfill, lca_declared_unit)
    df_emission.at[new_index, "endoflife_landfill"] = share_dross * func_lca_get_emission_factor("landfill", "dross")

    # recycling of metal materials
    # transport to recycling facility
    df_emission.at[new_index, "transport_recycler_metal"] = share_metal * func_lca_emissions_transport(lca_vehicle_to_recycler_metal, lca_distance_to_recycling_metal, lca_declared_unit)

    # sorting, shredding and cleaning of metals at recycling plant
    df_emission.at[new_index, "processing_metal_sorting"] = share_metal * func_lca_emissions_process(category="electricity", amount=lca_electricity_use_sorting_metal, origen=lca_emission_origen, source="mix")
    df_emission.at[new_index, "processing_metal_shredding"] = share_metal * func_lca_emissions_process(category="electricity", amount=lca_electricity_use_shredding_metal, origen=lca_emission_origen, source="mix")
    df_emission.at[new_index, "processing_metal_cleaning"] = share_metal * sum([
        func_lca_emissions_process(category="electricity", amount=lca_electricity_use_cleaning_metal, origen=lca_emission_origen, source="mix"), 
        func_lca_emissions_process(category="heat", amount=lca_heat_use_cleaning_metal, origen=lca_emission_origen, source="mix"),
        func_lca_emissions_process(category="water", amount=lca_water_use_cleaning_metal, origen="EU", source="groundwater"),
        func_lca_emissions_process(category="wastewater", amount=lca_wastewater_use_cleaning_metal, origen=lca_emission_origen, source="default"),
        func_lca_emissions_process(category="ressource", amount=lca_solvent_use_cleaning_metal, origen=lca_emission_origen, source="cleaning_agent")
    ])

    # melting down metals (depending on metal)
    df_emission.at[new_index, "processing_metal_melting"] = share_metal * func_lca_emissions_process_melting(df=df_result)

    # advantage due to secondary materials and energy generation
    #  Secondary Material (metal)
    df_emission.at[new_index, "avoided_electricity"] = lca_substitution_factor_electricity * -func_lca_emissions_process(category="electricity", amount=energy_incineration[0], origen=lca_emission_origen, source="mix")
    df_emission.at[new_index, "avoided_heat"] = lca_substitution_factor_heat * -func_lca_emissions_process(category="heat", amount=energy_incineration[1], origen=lca_emission_origen, source="mix")
    df_emission.at[new_index, "avoided_metal"]  = lca_substitution_factor_metal * -func_lca_emissions_avoided_material(df=df_result, relevant_category="metal")

    # Group emissions by category
    df_emission.at[new_index, "processing_metal_total"] = func_sum_columns_with_conditions(df=df_emission, row_label=new_index , include_substring="processing_metal", exclude_substring="total")
    df_emission.at[new_index, "endoflife_total"] = func_sum_columns_with_conditions(df=df_emission, row_label=new_index , include_substring="endoflife", exclude_substring="total")
    df_emission.at[new_index, "transport_total"] = func_sum_columns_with_conditions(df=df_emission, row_label=new_index , include_substring="transport", exclude_substring="total")
    df_emission.at[new_index, "avoided_total"] = func_sum_columns_with_conditions(df=df_emission, row_label=new_index , include_substring="avoided", exclude_substring="total")
    df_emission.at[new_index, "emissions_overall_per_declared_unit"] = func_sum_columns_with_conditions(df=df_emission, row_label=new_index , include_substring="total", exclude_substring="-")

    # Total emission per declared unit, by total waste weight (now and per year)
    material_weight_per_year = func_lca_get_weigth_per_year(material_weight_regular, material_frequency)
    
    df_emission.at[new_index, "emissions_overall_per_weight"] = waste_weigth * df_emission.at[new_index, "emissions_overall_per_declared_unit"]
    df_emission.at[new_index, "emissions_overall_per_year"] = material_weight_per_year * df_emission.at[new_index, "emissions_overall_per_declared_unit"]

    return df_emission

# function to calculate emissions for scenario: materials go to recycling
def func_lca_emissions_scenario_recycling(df: pd.DataFrame, df_emission: pd.DataFrame) -> pd.DataFrame: 
    
    # determine the new index to add the scenario as new row to the emission dataframe
    new_index = len(df_emission)

    # transport to recycling plant
    df_emission.at[new_index, "transport_recycler_plastic"] = func_lca_emissions_transport(lca_vehicle_to_recycler_plastic, lca_distance_to_recycler_plastic, lca_declared_unit)

    # recycling of materials
    # sorting, shredding, cleaning, drying and regranulation
    df_emission.at[new_index, "processing_plastic_sorting"] = func_lca_emissions_process(category="electricity", amount=lca_electricity_use_sorting_mixed, origen=lca_emission_origen, source="mix")

    share_metal, _ = func_lca_get_category_share(df_result, "metal") # share for remaining metal compared to 1 kg material
    share_plastic_init, _ = func_lca_get_category_share(df_result, "plastic") # share for remaining plastic compared to 1 kg material
    share_nonrecycable_init = 0 # TODO Abhängig machen zusammen mit share_plastic von den Kunststofffraktionen im Abfall, die nicht recycelbar sind (mit recycling Advancement aus material list, ggf. als weitere Spalte in df_result etablieren)

    # recycling of plastics (shredding, cleaning, drying, regranulation)
    df_emission.at[new_index, "processing_plastic_shredding"] = share_plastic_init * func_lca_emissions_process(category="electricity", amount=lca_electricity_use_shredding_plastic, origen=lca_emission_origen, source="mix")
    share_plastic = share_plastic_init * lca_share_io_shredding
    
    df_emission.at[new_index, "processing_plastic_cleaning"] = share_plastic * sum([
            func_lca_emissions_process(category="electricity", amount=lca_electricity_use_cleaning_plastic, origen=lca_emission_origen, source="mix"), 
            func_lca_emissions_process(category="heat", amount=lca_heat_use_cleaning_plastic, origen=lca_emission_origen, source="mix"),
            func_lca_emissions_process(category="water", amount=lca_water_use_cleaning_plastic, origen="EU", source="groundwater"),
            func_lca_emissions_process(category="wastewater", amount=lca_wastewater_use_cleaning_plastic, origen=lca_emission_origen, source="default"),
            func_lca_emissions_process(category="ressource", amount=lca_solvent_use_cleaning_plastic, origen=lca_emission_origen, source="cleaning_agent")
        ])    
    share_plastic *= lca_share_io_cleaning
    
    # depending on plastic type
    df_emission.at[new_index, "processing_plastic_drying"] = share_plastic * func_lca_emissions_processes_plastic(process="drying", df=df_result)
    df_emission.at[new_index, "processing_plastic_regranulation"] = share_plastic * func_lca_emissions_processes_plastic(process="regranulation", df=df_result)
    share_plastic *= lca_share_io_regranulation
    share_nonrecycable = share_nonrecycable_init + share_plastic_init * (1 - lca_share_io_overall)

    # recycling of metal materials
    # transport to recycling facility
    df_emission.at[new_index, "transport_recycler_metal"] = share_metal * func_lca_emissions_transport(lca_vehicle_to_recycler_metal, lca_distance_to_recycling_metal, lca_declared_unit)

    # sorting, shredding and cleaning of metals at recycling plant
    df_emission.at[new_index, "processing_metal_sorting"] = share_metal * func_lca_emissions_process(category="electricity", amount=lca_electricity_use_sorting_metal, origen=lca_emission_origen, source="mix")
    df_emission.at[new_index, "processing_metal_shredding"] = share_metal * func_lca_emissions_process(category="electricity", amount=lca_electricity_use_shredding_metal, origen=lca_emission_origen, source="mix")
    df_emission.at[new_index, "processing_metal_cleaning"] = share_metal * sum([
        func_lca_emissions_process(category="electricity", amount=lca_electricity_use_cleaning_metal, origen=lca_emission_origen, source="mix"), 
        func_lca_emissions_process(category="heat", amount=lca_heat_use_cleaning_metal, origen=lca_emission_origen, source="mix"),
        func_lca_emissions_process(category="water", amount=lca_water_use_cleaning_metal, origen="EU", source="groundwater"),
        func_lca_emissions_process(category="wastewater", amount=lca_wastewater_use_cleaning_metal, origen=lca_emission_origen, source="default"),
        func_lca_emissions_process(category="ressource", amount=lca_solvent_use_cleaning_metal, origen=lca_emission_origen, source="cleaning_agent")
    ])

    # melting down metals (depending on metal)
    df_emission.at[new_index, "processing_metal_melting"] = share_metal * func_lca_emissions_process_melting(df=df_result)

    # incineration of non-recycable materials
    # transport to waste-to-energy(wte)-plant
    df_emission.at[new_index, "transport_wte"] = share_nonrecycable * func_lca_emissions_transport(lca_vehicle_to_wte, lca_distance_to_wte, lca_declared_unit)

    # incineration of waste in waste-to-energy(wte)-plant
    df_emission.at[new_index, "endoflife_incineration"] = share_nonrecycable * func_lca_emissions_incineration(df_result, lca_declared_unit)
    energy_incineration = [x * share_nonrecycable for x in func_lca_energy_incineration(df_result, lca_declared_unit)]
    share_dross = lca_share_dross * share_nonrecycable # share for the remaining dross compared to 1 kg waste

    # transport of dross to landfill + enmissions from landfill
    df_emission.at[new_index, "transport_landfill"] = share_dross * func_lca_emissions_transport(lca_vehicle_to_landfill, lca_distance_to_landfill, lca_declared_unit)
    df_emission.at[new_index, "endoflife_landfill"] = share_dross * func_lca_get_emission_factor("landfill", "dross")

    # advantage due to secondary materials and energy generation
    df_emission.at[new_index, "avoided_electricity"] = lca_substitution_factor_electricity * -func_lca_emissions_process(category="electricity", amount=energy_incineration[0], origen=lca_emission_origen, source="mix")
    df_emission.at[new_index, "avoided_heat"] = lca_substitution_factor_heat * -func_lca_emissions_process(category="heat", amount=energy_incineration[1], origen=lca_emission_origen, source="mix")
    df_emission.at[new_index, "avoided_plastic"]  = lca_substitution_factor_plastic * -func_lca_emissions_avoided_material(df=df_result, relevant_category="plastic")
    df_emission.at[new_index, "avoided_metal"]  = lca_substitution_factor_metal * -func_lca_emissions_avoided_material(df=df_result, relevant_category="metal")

    # Group emissions by category
    df_emission.at[new_index, "processing_plastic_total"] = func_sum_columns_with_conditions(df=df_emission, row_label=new_index , include_substring="processing_plastic", exclude_substring="total")
    df_emission.at[new_index, "processing_metal_total"] = func_sum_columns_with_conditions(df=df_emission, row_label=new_index , include_substring="processing_metal", exclude_substring="total")
    df_emission.at[new_index, "endoflife_total"] = func_sum_columns_with_conditions(df=df_emission, row_label=new_index , include_substring="endoflife", exclude_substring="total")
    df_emission.at[new_index, "transport_total"] = func_sum_columns_with_conditions(df=df_emission, row_label=new_index , include_substring="transport", exclude_substring="total")
    df_emission.at[new_index, "avoided_total"] = func_sum_columns_with_conditions(df=df_emission, row_label=new_index , include_substring="avoided", exclude_substring="total")
    df_emission.at[new_index, "emissions_overall_per_declared_unit"] = func_sum_columns_with_conditions(df=df_emission, row_label=new_index , include_substring="total", exclude_substring="-")

    # Total emission per declared unit, by total waste weight (now and per year)
    material_weight_per_year = func_lca_get_weigth_per_year(material_weight_regular, material_frequency)
    
    df_emission.at[new_index, "emissions_overall_per_weight"] = waste_weigth * df_emission.at[new_index, "emissions_overall_per_declared_unit"]
    df_emission.at[new_index, "emissions_overall_per_year"] = material_weight_per_year * df_emission.at[new_index, "emissions_overall_per_declared_unit"]

    return df_emission

# compare emission values
def func_lca_compare_emissions(val_recycling: float, val_wte: float) -> tuple:
    """
    Compares emissions between recycling and waste incineration (WTE) and returns 
    both the difference and a sentence describing the comparison.

    Parameters:
        val_recycling (float): Emissions for recycling in kg.
        val_wte (float): Emissions for waste incineration in kg.

    Returns:
        tuple: A tuple containing:
            - float: The calculated difference in emissions (recycling perspective).
            - str: A sentence describing the comparison.
    """
    # Calculate the difference
    difference = val_wte - val_recycling
    
    # Create a sentence based on the difference
    if difference > 0:
        sentence = f"Die Emissionen für das Recycling sind {difference:.2f} kg niedriger als für die Abfallverbrennung."
    elif difference < 0:
        sentence = f"Die Emissionen für das Recycling sind {-difference:.2f} kg höher als für die Abfallverbrennung."
    else:
        sentence = "Die Emissionen für das Recycling und die Abfallverbrennung sind identisch."
    
    return difference, sentence

## Main script
# Step 1: Check hazourdous/non-hazourdous status via REACH-conformity
# 1.1 If waste is not conform to REACH, categorize as "Sondermüll"
if reach_conformity == "Nein":
    wertstoffscore = -10.0
    ws_category, ws_sentence = func_evaluateWS(wertstoffscore)
   
# Step 2: Check if fractions are for the most part potentially recycable
# 2.1 Filter materials dataframe for recycling advancement
df_materials_recyAdvance = df_materials[df_materials.recycling_advance > 0] #filter dataframe for all materials with the potential for recycling >0

# 2.2 Adding up all
for k in range(material_fraction_number):
  # Check if the kth entry of material_type is in a list of materials of dataframe
    if material_type[k] in df_materials_recyAdvance.abbreviation.tolist():
        # Add the kth value of material_share to percent_valuable_substance
        percent_valuable_substance += material_share[k]/100

if percent_valuable_substance < set_threshold_assessibility:
    wertstoffscore = 0
    ws_category, ws_sentence = func_evaluateWS(wertstoffscore)

# Step 3: Check if possible to sort fractions by different methods
# Initialize dataframes for output scores (1. for sorting, 2. for final results)
df_result = func_initializeResultDf(amount=material_fraction_number, label_columns=columns_materials)

# if waste consists only of one fraction, no sorting is necessary
if material_fraction_number == 1: 
    wertstoffscore = 89 #Einteilung als Recyling, werkstofflich, sortenrein
    ws_category, ws_sentence = func_evaluateWS(wertstoffscore)

    # result_sorting in Dataframe = 1
    df_result.at[0, label_material_result_sorting] = 1

# if waste consists of more than 1 fraction, sorting is necessary
elif material_fraction_number > 1 and material_fraction_number <= 5:
  
    # Store dataframes in list to use in loop
    list_df_sort = [df_sort_ferromagnetic, df_sort_eddycurrent, df_sort_density, df_sort_electrostatic, df_sort_sensorbased]
    list_sort_weight = [wt_ferromagnetic, wt_eddycurrent, wt_density, wt_electrostatic, wt_sensorbased]

    # Initialize dataframes for output scores (1. for sorting)
    df_result_sorting = func_initializeOutputDf(amount = material_fraction_number, label_columns=columns_sorting)

    # loop through each sorting method and obtain values for all material pairing
    for k in range(len(list_df_sort)):

        # Obtain values from function
        values_sort = func_checkSorting(amount=material_fraction_number, df=list_df_sort[k])

        # Add values for each sorting method to result dataframe and change index names
        df_result_sorting[columns_sorting[k+1]] = values_sort

    # Check if one material can be sorted completly from any other (= 1, sortenrein)
    for k in range(material_fraction_number):

        # Get material string 
        material_name = material_type[k]

        # Check which row contains material name
        matching_row_indices = func_get_indices_for_material(df=df_result_sorting, column=label_material_pairing, material=material_name)
        st.write(matching_row_indices)

        # Check if all relevant rows for a material contain at least a 1
        value_to_find = 1
        list_check_clear_sort = []

        for l in range(len(matching_row_indices)):
            contains_value = value_to_find in df_result_sorting.iloc[matching_row_indices[l], :].values
            list_check_clear_sort.append(contains_value)

        # if the material can be sorted from other materials: 
        if all(list_check_clear_sort) == True:
            # res_sort = 1
            df_result.at[k, label_material_result_sorting] = 1

        # if the material can NOT be sorted completly from other materials:  
        elif all(list_check_clear_sort) == False:
            # res_sort = weigthed average of all results by sorting method for specific material
            for l in range(len(list_check_clear_sort)):
                
                # get material pairing which cannot be sorted completly (!=1)
                if list_check_clear_sort[l] == False:
                    # extract values from dataframe column
                    values_to_average = df_result_sorting.iloc[:, matching_row_indices[l]].tolist()
                    # multiply each value with the corresponding weight
                    weighted_values_to_average = [a * b for a, b in zip(values_to_average, list_sort_weight)]
                    # calculate mean and store in result dataframe
                    df_result.at[k, label_material_result_sorting] = sum(weighted_values_to_average) / len(weighted_values_to_average)

    wertstoffscore = 74 #Einteilung als Recyling, werkstofflich, gemischt
    ws_category, ws_sentence = func_evaluateWS(wertstoffscore)

st.write(df_result_sorting)
st.write(df_result)

## get location data from Recycler #TODO WeSort: Hier den Algorithmus zur Verknüpfung mit dem Recycler einfügen. 
# Habe bereits die ermittelten Daten zu Längen- und Breitengrad des Abfallursprungs als list bereitgestellt. 
# Als Output sollte u.a. die Transportdistanz vom Unternehmen zum Recycler (lca_distance_to_recycler) angegeben werden. Diese wird unten für die life cycle analysis benötigt.
#company_coordinates = st.session_state.coordinates_data

lca_distance_to_recycler_plastic = 10 #TODO WeSort: Hier mit dem Ausgabewert für die Distance ersetzen.
contact_recycler_name = "Multiport GmbH" # TODO WeSort: Hier Beispiel austauschen mit gewähltem Recycler
contact_recycler_address = "Ernst-Grube-Straße 1, 06406 Bernburg" # TODO WeSort: Hier Beispiel austauschen mit gewähltem Recycler
contact_recycler_phone = "+49 3471 64040" # TODO WeSort: Hier Beispiel austauschen mit gewähltem Recycler
contact_recycler_mail = "de-ves-multiport-bernburg@veolia.com" # TODO WeSort: Hier Beispiel austauschen mit gewähltem Recycler

## life cycle analysis (lca) and comparison of current and proposed waste treatment

# Initialize dataframe for calculated emissions
df_emission = pd.DataFrame(columns = columns_emissions)

# Scenario 1: Calulation of lca for current waste treatment (assumption waste-to-energy)
# Use function to calculate emissions in this scenario
df_emission = func_lca_emissions_scenario_wte(df=df_result, df_emission=df_emission)

# Scenario 2. Calculation of lca of future waste treatment (recycling)
df_emission = func_lca_emissions_scenario_recycling(df=df_result, df_emission=df_emission)

# Compare emissions from both scenario
emission_recyling_per_weight = df_emission.at[1, "emissions_overall_per_weight"] # emission in t CO2 / t material
emission_wte_per_weight = df_emission.at[0, "emissions_overall_per_weight"] # emission in t CO2 / t material
emission_comparison_per_weight, emission_sentence_per_weight = func_lca_compare_emissions(df_emission.at[1, "emissions_overall_per_weight"], df_emission.at[0, "emissions_overall_per_weight"])
emission_comparison_per_year, emission_sentence_per_year = func_lca_compare_emissions(df_emission.at[1, "emissions_overall_per_year"], df_emission.at[0, "emissions_overall_per_year"])

# Title Section
st.title("PlastIQ")
st.subheader("Dein Wertstoff-Score (WS)")

# Main Content
col1, col2 = st.columns([1, 2])  # Adjust proportions to balance text and map

# Left column: Score and recommendation
with col1:
    st.metric(label="WS", value=f"{wertstoffscore} %", help="Der Wertstoff-Score wurde auf Basis deiner Eingaben ermittelt. Er gibt an, welche Verwertungsmethode dafür geeignet erscheint.")
    st.write(f"Dein Wertstoff erreicht den Score {ws_category}. {ws_sentence}")

    # Recycler Recommendation
    st.write("### Wir empfehlen dir folgenden Recycler:")
    st.write(f"**{contact_recycler_name}**")
    st.write(contact_recycler_address)
    st.write(f"**Kontakt:** {contact_recycler_phone}")
    st.write(f"[{contact_recycler_mail}](mailto:{contact_recycler_mail})")

# Right column: Map (placeholder)
with col2:
    # Replace this with your dynamic map if needed
    st.map()  # Placeholder for location/map integration TODO WeSort: Hier die Karte mit den Standorten des Unternehmens, des Recyclers und der Fahrtstrecke anzeigen lassen

# Economic and Ecological Comparison
st.divider()  # Horizontal divider

col3, col4 = st.columns([3,1])

# Ecological Section
with col3:
    st.header("Ökologisch")
    st.write(f"THG-Emission Recycling (erwartet): **{emission_recyling_per_weight.round(1)} t CO₂e/t**")
    st.write(f"THG-Emission Entsorgung bisher (geschätzt): **{emission_wte_per_weight.round(1)} t CO₂e/t**")
    st.write(f"Deine Wertstoffmenge: **{waste_weigth} t**")
    st.subheader(f"Einsparung: **{emission_comparison_per_weight.round(1)} t CO₂e**")

# Economical Section
#with col4:
#    st.header("Ökonomisch")
#    st.write("Erlöse für Recycling (erwartet): **500 €/t**")
#    st.write("Kosten für Entsorgung bisher (geschätzt): **-100 €/t**")
#    st.write("Ihre Wertstoffmenge: **2,5 t**")
#    st.subheader("Vorteil: **1.500 €**")