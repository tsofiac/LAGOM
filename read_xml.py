import xml.etree.ElementTree as ET
import csv

def xml_to_csv(xml_file, csv_file):
    # Parse the XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Open a file for writing
    with open(csv_file, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)

        # Write the header row based on XML tags
        # Assuming all XML entries have the same structure
        first_record = root[0]
        header = [elem.tag for elem in first_record]
        writer.writerow(header)

        # Write data rows
        for record in root:
            row = [elem.text for elem in record]
            writer.writerow(row)

# Example usage
xml_to_csv('dataset/raw_data/drugbank_full_database.xml', 'drugbank_full.csv')