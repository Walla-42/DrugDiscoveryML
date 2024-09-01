#!/bin/bash
cd ..
wget http://www.yapcwsoft.com/dd/padeldescriptor/PaDEL-Descriptor.zip


unzip PaDEL-Descriptor.zip -d ./PaDEL-Descriptor

rm PaDEL-Descriptor.zip

echo "PaDEL-Descriptor downloaded and extracted to ./PaDEL-Descriptor"

# Read the output directory from the file
output_dir=$(cat ./PaDEL-Descriptor/output_dir.txt)
output_file="$output_dir/descriptors_output.csv"


# Check if the directory exists
if [ ! -d "$output_dir" ]; then
  echo "Directory $output_dir does not exist. Please check the output_dir.txt file."
  exit 1
fi

# Execute the Java program
java -Xms1G -Xmx1G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar \
 -removesalt -standardizenitro -fingerprints \
 -descriptortypes ./PaDEL-Descriptor/descriptors.xml \
 -dir "$output_dir" -file "$output_file"




