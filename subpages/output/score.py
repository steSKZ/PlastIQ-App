## Import libraries
import streamlit as st
import numpy as np
import pandas as pd
from itertools import chain

## Define variables
# Initialize variables
wertstoffscore = 0 # variable for the final score of the waste material
percent_valuable_substance = 0.0 #set up variable for the amount of valuable substance in the waste
waste_type = [] # list of waste types added during the input
waste_share = [] # list of corresponding waste percentages added during the input

# Parameters changing model
set_threshold_assessibility = 0.8 # threshold for Step 2 regarding the assessibility of the recycling fraction
wt_ferromagnetic = 1 # weigth of sorting method "ferromagnetic", between 0 and 1
wt_eddycurrent = 1 # weigth of sorting method "eddycurrent", between 0 and 1
wt_density = 1 # weigth of sorting method "density", between 0 and 1
wt_electrostatic = 1 # weigth of sorting method "electrostatic", between 0 and 1

# Default LCA data
lca_declared_unit = 1 # Declared unit, usually 1 kg

lca_substitution_factor_plastics = 0.5 # substitution factor for avoided emissions for plastics # TODO 
lca_substitution_factor_metal = 0.8 # substitution factor for avoided emissions for metal # TODO 
lca_substitution_factor_electricity = 1 # substitution factor for avoided emissions for electricity # TODO
lca_substitution_factor_heat = 1 # substitution factor for avoided emissions for heat # TODO 

lca_distance_to_wte = 10 # Distance to waste-to-energy-plant in km TODO
lca_distance_to_landfill = 10 # Distance from wte-plant to landfill in km TODO
lca_distance_to_recycling_metal = 10 # Distance from wte-plant to recycling facility for metal in km TODO

lca_efficiency_wte_electric = 0.113
lca_efficiency_wte_heat = 0.33

lca_electricity_use_sorting_metal = 0.3 # electricity use for sorting metal per kg, in kWh TODO
lca_electricity_use_shredding_metal = 0.3 # electricity use for shredding metal per kg, in kWh TODO
lca_electricity_use_cleaning_metal = 0.3 # electricity use for cleaning metal per kg, in kWh TODO
lca_electricity_use_melting_metal = 1 # electricity use for melting metal per kg, in kWh TODO
lca_heat_use_cleaning_metal = 1 # heat energy use for cleaning metal per kg, in kWh TODO
lca_water_use_cleaning_metal = 1 # water use for cleaning metal per kg in kg Water TODO
lca_wastewater_use_cleaning_metal = 1 # wastewater for cleaning metal per kg in kg wastewater TODO
lca_solvent_use_cleaning_metal = 1 # solvent use for clening metal per kg in kg solvent TODO

lca_energy_use_sorting_mixed = 0.3 # electricity use for sorting mixed materials per kg in kWh TODO

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
columns_sorting = []
# for dataframe emissions
columns_emissions = [
    "processing_total",
    "processing_plastic_total", "processing_plastic_sorting", "processing_plastic_shredding", "processing_plastic_cleaning", "processing_plastic_drying", "processing_plastic_regranulation", 
    "processing_metal_total", "processing_metal_sorting", "processing_metal_shredding", "processing_metal_cleaning", "processing_metal_melting",
    "endoflife_total", "endoflife_incineration", "endoflife_landfill",
    "transport_total", "transport_recycler_plastic", "transport_recycler_metal", "transport_wte", "transport_landfill",
    "avoided_total", "avoided_electricity", "avoided_heat", "avoided_plastic", "avoided_metal"
]   

## Define functions - general
# Function to locate WS to a class
def func_evaluateWS(wertstoffscore: float):
    ws_rounded = round(wertstoffscore, 1)
    step = 0.1
    if -10 <= ws_rounded < 0: #Klasse G: Sondermüll
      ws_category = "H"
    elif 0 <= ws_rounded < 50: #Klasse F: Energetische Verwertung
      ws_category = "F"
    elif 50 <= ws_rounded < 60: #Klasse E: Recycling, chemisch
      ws_category = "E"
    elif 60 <= ws_rounded < 75: #Klasse D: Recycling, werkstofflich, gemischt
      ws_category = "D"
    elif 75 <= ws_rounded < 90: #Klasse C: Recyling, werkstofflich, sortenrein
      ws_category = "C"
    elif 90 <= ws_rounded < 95: #Klasse B: Wiederverwendung, niederwertig
      ws_category = "B"
    elif 95 <= ws_rounded <= 100: #Klasse A: Wiederverwendung, gleichwertig
      ws_category = "A"
    else: #Fehlerhafte Eingabe
      ws_category = "Eingabe konnte nicht verarbeitet werden."
    return ws_category


# Initialize dataframe for output scores
def func_initializeOutputDf(count, waste_type):
    # Loop through each entry and compare it with every other entry, while avoiding duplicates
    df_new = pd.DataFrame()
    
    # Add column for each material pairing
    for k in range(count):
        for l in range(k+1, count):
            
            # Specify the row and column indices
            first_material = waste_type[k]
            second_material = waste_type[l]

            # Define material pairing
            pair_string = f"{first_material}_{second_material}"

            # Add to dataframe
            df_new[pair_string] = []

    # Output new dataframe
    return df_new

# Initialize dataframe for result 
def func_initializeResultDf(label_columns: list, type: list, share: list, df_materials: pd.DataFrame):
    
    # New dataframe with material as column name
    df_new = pd.DataFrame(columns= label_columns)
    # Add first column "material type"
    df_new[label_material_type] = type
    # Lookup category of material
    material_category = []

    for k in range(len(type)):
        
        # lookup category from material dataframe
        category = df_materials.loc[df_materials["abbreviation"].str.fullmatch(type[k], case=False, na=False), "category"]
        material_category.append(category.iloc[0])

    # Add second column with material category
    df_new[label_material_category] = material_category

    # Add third column with material share
    df_new[label_material_share] = share            

    return df_new

# Function to check all available materials for their sorting options
def func_checkSorting(count, df_name):
    
    # Initialize list with all result values
    values_sort = []

    # Loop through each entry and compare it with every other entry, while avoiding duplicates
    for k in range(count):
        for l in range(k+1, count):
              
          # Specify the row and column indices
          row_label = waste_type[k]
          column_label = waste_type[l]

          # Get and return value at indices from matrix
          result = df_name.loc[row_label, column_label]

          # Add to list
          values_sort.append(float(result))
    
    # Output list with result values
    return values_sort

## Define functions - LCA
# Get a specific emission factor from background_data
def func_lca_get_emission_factor(df: pd.DataFrame, label_category: str, label_use: str):
    # Filter df by specific category
    df_filtered = df[df["category"] == label_category]

    # lookup emission factor specific to the label_use
    emission_factor = float(df_filtered.loc[df_filtered["use for"].str.contains(label_use, case=False, na=False), "GWP100"])

    return emission_factor

# Get emissions from transport
def func_lca_emissions_transport(df: pd.DataFrame, label_use: str, distance: float, payload: float):
    # get emission factor
    emission_factor = func_lca_get_emission_factor(df, label_category="transport", label_use=label_use)
    # calculate emissions with payload (in ton) * distance (in km) * emission factor (in kg CO2e/(t*km))
    emission_transport = payload * distance * emission_factor
    
    return emission_transport

# get emissions from waste incineration
def func_lca_emissions_incineration(fractions, df_emissions: pd.DataFrame, df_result: pd.DataFrame, df_materials: pd.DataFrame, weight: float):
    
    # Initialize values for resulting emissions, electric and heat energy
    emissions_incineration = 0 #in kg CO2e/kg waste
    electric_energy_incineration = 0 #in kWh 
    heat_energy_incineration = 0 #in kWh

    # for every waste fraction
    for k in range(fractions):
        
        # check if material is plastic (has heating value and emission) or metal (nor heating value nor emission)
        material_category = df_result.at["category", df_result.columns[k]]
        
        if material_category == "plastic":
            # get waste name and percentage from df_result in a list [name, percentage]
            waste_data = [df_result.columns[k], df_result.at["share", df_result.columns[k]]]
            st.write(waste_data)
            # lookup name in df_emission, get emission_factor and add to list
            emission_factor = func_lca_get_emission_factor(df_emissions, "incineration", waste_data[0])

            # lookup heating value from df_material
            lower_heating_value = float(df_materials.loc[df_materials["abbreviation"].str.fullmatch(waste_data[0], case=False, na=False), "lower_heating_value_MJ_per_kg"])

            # calculate electric and heat energy [kWh] from heating value [J/kg], weight [kg], material_share [-] and efficiency [-]
            electric_energy_incineration += lower_heating_value * weight * waste_data[1]/100 * lca_efficiency_wte_electric / 3.6
            heat_energy_incineration += lower_heating_value * weight * waste_data[1]/100 * lca_efficiency_wte_heat / 3.6

            # add to exisiting emissions and energy
            emissions_incineration += weight * waste_data[1]/100 * emission_factor
            st.write(emissions_incineration)
            st.write(electric_energy_incineration)
    return [emissions_incineration, [electric_energy_incineration, heat_energy_incineration]]

# get emissions from process (categories: electricity, heat, ressources)
def func_lca_emissions_process(df: pd.DataFrame, category: str, amount: float, origen: str, source: str):
    # define label for emission factor search 
    label_use = f"{category}_{origen}_{source}"
    # look up emission factor
    emission_factor = func_lca_get_emission_factor(df, label_category = category, label_use = label_use)
    # calculate emission (kg CO2e/kWh) from electrical energy (in kWh) and emission factor
    emission_process= amount * emission_factor

    return emission_process

# get avoided emissions for secondary materials 1. plastic, 2. metal, depending on material in stream
def func_lca_emission_avoided_material(df_result: pd.DataFrame, df_emissions: pd.DataFrame, relevant_category: str):
    
    # loop through every type of material
    for k in df_result.shape[0]:

        # get material type, category and share from df_result
        material_name = df_result.columns.tolist[k]
        material_category = df_result.loc["category", material_name]
        material_share = df_result.loc["share (%)", material_name]
        
        # Summe der material-Anteile der entsprechnedne Kategorie
        category_share = float(df_result[df_result["category"] == relevant_category].loc["share (%)"].sum)

        if material_category == relevant_category:

            emission_avoided_material =  (material_share/category_share) * func_lca_emissions_process()



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

## Main script
# Step 1: Check hazourdous/non-hazourdous status via REACH-conformity
reach_conformity = st.session_state.key_dict_product_quality["input_wertstoff_reach"]

# 1.1 If waste is not conform to REACH, categorize as "Sondermüll"
if reach_conformity == "Nein":
    wertstoffscore = -10.0
    ws_category = func_evaluateWS(wertstoffscore)
   
# Step 2: Check if fractions are for the most part potentially recycable
# 2.1 Load background data and set up variables
df_materials = pd.read_excel(file_path_background, sheet_name = "list_material")
df_materials_recyAdvance = df_materials[df_materials.recycling_advance > 0] #filter dataframe for all materials with the potential for recycling >0
waste_fractions_number = st.session_state.key_dict_product["input_waste_fraction_number"]

for k in range(waste_fractions_number):
    waste_type.append(st.session_state.key_dict_product[f"input_wertstoff_typ_{k}"])
    waste_share.append(st.session_state.key_dict_product[f"input_wertstoff_anteil_{k}"])

# 2.2 Adding up all
for k in range(waste_fractions_number):
  # Check if the kth entry of waste_type is in a list of materials of dataframe
    if waste_type[k] in df_materials_recyAdvance.abbreviation.tolist():
        # Add the kth value of waste_share to percent_valuable_substance
        percent_valuable_substance += waste_share[k]/100

if percent_valuable_substance < set_threshold_assessibility:
    wertstoffscore = 0
    ws_category = func_evaluateWS(wertstoffscore)

# Step 3: Check if possible to sort fractions by different methods

# if waste consists only of one fraction, no sorting is necessary
if waste_fractions_number == 1: 
    wertstoffscore = 89 #Einteilung als Recyling, werkstofflich, sortenrein
    ws_category = func_evaluateWS(wertstoffscore)

# if waste consists of more than 1 fraction, sorting is necessary
elif waste_fractions_number > 1 and waste_fractions_number <= 5:
    
    # load relevant background data
    # Magnetabscheidung
    df_sort_ferromagnetic = pd.read_excel(file_path_background, sheet_name = "sort_ferromagnetic", index_col=0)
    # Wirbelstromsortierung
    df_sort_eddycurrent = pd.read_excel(file_path_background, sheet_name = "sort_eddycurrent", index_col=0)
    # Dichtesortierung
    df_sort_density = pd.read_excel(file_path_background, sheet_name = "sort_density", index_col=0)
    # Elektrostatische Sortierung
    df_sort_electrostatic = pd.read_excel(file_path_background, sheet_name = "sort_electrostatic", index_col=0)

    # Store dataframes in list to use in loop
    list_df_sort = [df_sort_ferromagnetic, df_sort_eddycurrent, df_sort_density, df_sort_electrostatic]
    list_df_sort_name = ["sort_ferromagnetic", "sort_eddycurrent", "sort_density", "sort_electrostatic"]
    list_sort_weight = [wt_ferromagnetic, wt_eddycurrent, wt_density, wt_electrostatic]

    # Initialize dataframes for output scores (1. for sorting, 2. for final results)
    df_result_sorting = func_initializeOutputDf(waste_fractions_number, waste_type)
    df_result = func_initializeResultDf(columns_materials, waste_type, waste_share, df_materials)

    # loop through each sorting method and obtain values for all material pairing
    for k in range(len(list_df_sort)):

        # Obtain values from function
        values_sort = func_checkSorting(waste_fractions_number,list_df_sort[k])

        # Add values for each sorting method to result dataframe and change index names
        df_result_sorting.loc[k] = values_sort
    
    # Change indexes of sorting method
    df_result_sorting.index = list_df_sort_name

    # add row for final sorting results to result dataframe
    #df_result[label_material_result_sorting] = [np.nan] * waste_fractions_number

    # Check if one material can be sorted completly from any other (= 1, sortenrein)
    for k in range(waste_fractions_number):

        # Get material string 
        material_name = waste_type[k]
        # Check which column contains material name
        matching_col_indices = [i for i, col in enumerate(df_result_sorting.columns) if material_name in col]
        # Check if all relevant columns for a material contain at least a 1
        value_to_find = 1
        list_check_clear_sort = []

        for l in range(len(matching_col_indices)):
            contains_value = value_to_find in df_result_sorting.iloc[:, matching_col_indices[l]].values
            list_check_clear_sort.append(contains_value)

        # if the material can be sorted from other materials: 
        if all(list_check_clear_sort) == True:
            # res_sort = 1
            df_result.at["res_sort", material_name] = 1

        # if the material can NOT be sorted completly from other materials:  
        elif all(list_check_clear_sort) == False:
            # res_sort = weigthed average of all results by sorting method for specific material
            for l in range(len(list_check_clear_sort)):
                
                # get material pairing which cannot be sorted completly (!=1)
                if list_check_clear_sort[l] == False:
                    # extract values from dataframe column
                    values_to_average = df_result_sorting.iloc[:, matching_col_indices[l]].tolist()
                    # multiply each value with the corresponding weight
                    weighted_values_to_average = [a * b for a, b in zip(values_to_average, list_sort_weight)]
                    # calculate mean and store in result dataframe
                    df_result.at["res_sort", material_name] = sum(weighted_values_to_average) / len(weighted_values_to_average)

st.write(df_result) #TODO

## get location data from Recycler #TODO WeSort: Hier den Algorithmus zur Verknüpfung mit dem Recycler einfügen. 
# Habe bereits die ermittelten Daten zu Längen- und Breitengrad des Abfallursprungs als list bereitgestellt. 
# Als Output sollte u.a. die Transportdistanz vom Unternehmen zum Recycler (lca_distance_to_recycler) angegeben werden. Diese wird unten für die life cycle analysis benötigt.
#company_coordinates = st.session_state.coordinates_data

lca_distance_to_recycler = 10 #TODO WeSort: Hier mit dem Ausgabewert für die Distance ersetzen.

## life cycle analysis (lca) and comparison of current and proposed waste treatment
# Get emission data
df_lca_emission_data = pd.read_excel(file_path_background, sheet_name = "lca_calculation")
# Get weight of materials
waste_weigth = st.session_state.key_dict_product_amount["input_wertstoff_menge"]
# Get weight of material per year
material_weight_regular = st.session_state.key_dict_product_amount["input_haeufigkeit_menge"]
material_frequency = st.session_state.key_dict_product_amount["input_haeufigkeit_turnus"]
material_weight_per_year = func_lca_get_weigth_per_year(material_weight_regular, material_frequency)

# 1. lca of current waste treatment (waste incineration) for 1kg of waste
# 1.1. Transport to waste-to-energy(wte)-plant
emission_current_transport_wte = func_lca_emissions_transport(df_lca_emission_data, "transport_lkw_22", lca_distance_to_wte, lca_declared_unit)

# 1.2. Incineration of waste in waste-to-energy(wte)-plant
[emission_current_incineration, energy_current_incineration] = func_lca_emissions_incineration(waste_fractions_number, df_lca_emission_data, df_result, df_materials, lca_declared_unit)

# 1.3. Sorting process after incineration between metals (recycling) and dross (landfill)

dross_share = 0.2 # share for the remaining dross compared to 1 kg waste # TODO
metal_share = 0.5 # share for remaining metal compared to 1 kg waste # TODO

# 1.4 Transport of dross to landfill + enmissions from landfill
emission_current_transport_landfill = dross_share * func_lca_emissions_transport(df_lca_emission_data, "transport_lkw_22", lca_distance_to_landfill, lca_declared_unit)
emission_current_dross_landfill = dross_share * func_lca_get_emission_factor(df_lca_emission_data, "landfill", "dross")

# 1.5. Recycling of metal materials
# 1.5.1. Transport to recycling facility
emission_current_transport_metal = metal_share * func_lca_emissions_transport(df_lca_emission_data, "transport_lkw_22", lca_distance_to_recycling_metal, lca_declared_unit)

# 1.5.2. Sorting, shredding and cleaning of metals at recycling plant
emsission_current_sorting_metal = metal_share * func_lca_emissions_process(df_lca_emission_data, "electricity", lca_electricity_use_sorting_metal, "DE", "mix")
emission_current_shredding_metal = metal_share * func_lca_emissions_process(df_lca_emission_data, "electricity", lca_electricity_use_shredding_metal, "DE", "mix")
emission_current_cleaning_metal = metal_share * sum([
    func_lca_emissions_process(df_lca_emission_data, "electricity", lca_electricity_use_cleaning_metal, "DE", "mix"), 
    func_lca_emissions_process(df_lca_emission_data, "heat", lca_heat_use_cleaning_metal, "DE", "mix"),
    func_lca_emissions_process(df_lca_emission_data, "water", lca_water_use_cleaning_metal, "EU", "groundwater"),
    func_lca_emissions_process(df_lca_emission_data, "wastewater", lca_wastewater_use_cleaning_metal, "DE", "default"),
    func_lca_emissions_process(df_lca_emission_data, "ressource", lca_solvent_use_cleaning_metal, "DE", "cleaning_agent")
])

# 1.5.3. Melting down metals (depending on metal)
emissions_current_melting_metal = metal_share * 0 # TODO

# 1.6. Advantage due to secondary materials and energy generation
#  Secondary Material (metal)
emission_current_avoided_electricity = lca_substitution_factor_electricity * -func_lca_emissions_process(df_lca_emission_data, "electricity", energy_current_incineration[0], "DE", "mix")
emission_current_avoided_heat = lca_substitution_factor_heat * -func_lca_emissions_process(df_lca_emission_data, "heat", energy_current_incineration[1], "DE", "mix")
emission_current_avoided_secondary_metal = lca_substitution_factor_metal * -func_lca_emissions_process(df_lca_emission_data, "production-metal", energy_current_incineration[1], "DE", "mix") # TODO depending on metal

# Group emissions by category
emission_cat_current_transport = [emission_current_transport_wte, emission_current_transport_landfill, emission_current_transport_metal]
emission_cat_current_endoflife = [emission_current_incineration, emission_current_dross_landfill]
emission_cat_current_process = [emsission_current_sorting_metal, emission_current_shredding_metal, emission_current_cleaning_metal, emissions_current_melting_metal]
emission_cat_current_avoided = [emission_current_avoided_electricity, emission_current_avoided_heat, emission_current_avoided_secondary_metal]

# Total emission per declared unit, by total waste weight (now and per year)
emission_current_per_declared_unit = sum(chain(emission_cat_current_transport, emission_cat_current_endoflife, emission_cat_current_process, emission_cat_current_avoided))
emission_current_per_weight = waste_weigth * emission_current_per_declared_unit
emission_current_per_year = material_weight_per_year * emission_current_per_declared_unit

## 2. lca of future waste treatment (recycling) for 1kg of waste
# 2.1. Transport to recycler (with specific location data) 
emission_future_transport_recycler = func_lca_emissions_transport(df_lca_emission_data, "transport_lkw_22", lca_distance_to_recycler, lca_declared_unit)

# 2.2. Recycling of materials
# 2.2.1. sorting, shredding, cleaning, drying and regranulation
emission_future_sorting = func_lca_emissions_process(df_lca_emission_data, "electricity", lca_energy_use_sorting_mixed, "DE", "mix")

share_material_plastics = 0.8 # TODO
share_material_metal = 0.2 # TODO
share_material_nonrecycable = 0.1

# 2.2.2 Recycling of plastics (shredding, cleaning, drying, regranulation)
emission_future_shredding_plastic = share_material_plastics * 0 # TODO
emission_future_cleaning_plastic = share_material_plastics * 0 # TODO
# depending on plastic type
emission_future_drying_plastic = share_material_plastics * 0 # TODO
emissin_future_regranulation_plastic = share_material_plastics * 0 # TODO

# 2.2.3. Recycling of metal (transport, shredding, cleaning, melting)
emission_future_transport_metal = share_material_metal * func_lca_emissions_transport(df_lca_emission_data, "transport_lkw_22", lca_distance_to_recycling_metal, lca_declared_unit)

emission_future_shredding_metal = share_material_metal * func_lca_emissions_process(df_lca_emission_data, "electricity", lca_electricity_use_shredding_metal, "DE", "mix")
emission_furture_cleaning_metal = share_material_metal * sum([
    func_lca_emissions_process(df_lca_emission_data, "electricity", lca_electricity_use_shredding_metal, "DE", "mix"), 
    func_lca_emissions_process(df_lca_emission_data, "heat", lca_electricity_use_shredding_metal, "DE", "mix"),
    func_lca_emissions_process(df_lca_emission_data, "water", lca_electricity_use_shredding_metal, "DE", "mix"),
    func_lca_emissions_process(df_lca_emission_data, "wastewater", lca_electricity_use_shredding_metal, "DE", "mix"),
    func_lca_emissions_process(df_lca_emission_data, "ressource", lca_electricity_use_shredding_metal, "DE", "mix")
])
emissions_future_melting_metal = share_material_metal * 0 # TODO

# 2.3. Incineration of non-recycable materials
# 2.3.1. Transport of non-recycable materials to waste-to-energy-plant
emission_future_transport_wte = share_material_nonrecycable * func_lca_emissions_transport(df_lca_emission_data, "transport_lkw_22", lca_distance_to_wte, lca_declared_unit)

# 2.3.2. Incineration of non-recycable plastics in waste-to-energy-plant
[emission_future_incineration, energy_future_incineration] = share_material_nonrecycable * func_lca_emissions_incineration(waste_fractions_number, df_lca_emission_data, df_result, df_materials, lca_declared_unit) # TODO change df_results to only include plastic materials (non-recycables)
st.write([emission_future_incineration, energy_future_incineration])

dross_share_future = share_material_nonrecycable * dross_share # share for the remaining dross compared to 1 kg waste

# 2.3.3 Transport of dross to landfill + enmissions from landfill
emission_future_transport_landfill = dross_share_future * func_lca_emissions_transport(df_lca_emission_data, "transport_lkw_22", lca_distance_to_landfill, lca_declared_unit)
emission_future_dross_landfill = dross_share_future * func_lca_get_emission_factor(df_lca_emission_data, "landfill", "dross")

# 2.4. Advantage due to secondary materials and energy generation (negative emissions)
emission_future_avoided_electricity = lca_substitution_factor_electricity * -func_lca_emissions_process(df_lca_emission_data, "electricity", energy_future_incineration[0], "DE", "mix")
emission_future_avoided_heat = lca_substitution_factor_heat * -func_lca_emissions_process(df_lca_emission_data, "heat", energy_future_incineration[1], "DE", "mix")
emission_future_avoided_secondary_plastic = 0 # TODO
emission_future_avoided_secondary_metal = 0 # TODO

# Group emissions by category
emission_cat_future_transport = [emission_future_transport_recycler, emission_future_transport_metal, emission_future_transport_wte, emission_future_transport_landfill]
emission_cat_future_endoflife= [emission_future_incineration, emission_future_dross_landfill]
emission_cat_future_process_plastic = [share_material_plastics * emission_future_sorting, emission_future_shredding_plastic, emission_future_cleaning_plastic, emission_future_drying_plastic, emissin_future_regranulation_plastic]
emission_cat_future_process_metal = [share_material_metal * emission_future_sorting, emission_future_shredding_metal, emission_furture_cleaning_metal, emissions_future_melting_metal]
emission_cat_future_process = emission_cat_future_process_plastic + emission_cat_future_process_metal
emission_cat_future_avoided = [emission_future_avoided_electricity, emission_future_avoided_heat, emission_future_avoided_secondary_plastic, emission_future_avoided_secondary_metal]

# Total emission per declared unit, by total waste weight (now and per year)
emission_future_per_declared_unit = sum(chain(emission_cat_future_transport, emission_cat_future_endoflife, emission_cat_future_process, emission_cat_future_avoided))
emission_future_per_weight = waste_weigth * emission_future_per_declared_unit
emission_future_per_year = material_weight_per_year * emission_future_per_declared_unit

st.write(emission_current_per_year)
st.write(emission_future_per_year)
