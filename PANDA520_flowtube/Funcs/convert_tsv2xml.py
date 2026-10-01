import csv
import xml.etree.ElementTree as ET
import os

def tsv_to_custom_xml(tsv_file, output_folder, xml_filename):
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)
    
    # Full path to the XML file
    xml_file = os.path.join(output_folder, xml_filename)
    
    # Create the root element with the namespace
    root = ET.Element("mechanism", xmlns="https://mcm.york.ac.uk/MCM")
    
    # Create the species_defs element
    species_defs = ET.SubElement(root, "species_defs")
    
    # Read the TSV file
    with open(tsv_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')  # Read as dictionary using the header
        species_number = 1  # Initialize species number counter
        for row in reader:
            if not row.get("Name") or not row.get("Smiles"):
                continue
            # Create a species element
            species = ET.SubElement(species_defs, "species", {
                "species_number": f"s{species_number}",
                "species_name": row.get("Name", "Unknown")  # Use 'name' column
            })
            
            # Add SMILES element
            smiles = ET.SubElement(species, "smiles")
            smiles.text = row.get("Smiles", "Unknown")  # Use 'smiles' column
            
            species_number += 1  # Increment species number

    tree = ET.ElementTree(root)
    with open(xml_file, "w", encoding="utf-8") as file:
        file.write('<?xml version="1.0" encoding="UTF-8"?>\n')  # Add XML declaration manually
        file.write(
            ET.tostring(
                root, encoding="unicode", method="xml"
            ).replace("><", ">\n<")  # Ensure line breaks between elements
        )
    

