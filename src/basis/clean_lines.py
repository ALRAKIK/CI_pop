
# Python code to remove trailing spaces from each line in the file

file_path = './tmp/Basis_normalized'

# Read the file, process each line, and write back to the same file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Remove trailing spaces from each line
cleaned_lines = [line.rstrip() + '\n' for line in lines]

# Write the cleaned lines back to the file
with open(file_path, 'w') as file:
    file.writelines(cleaned_lines)

# Confirming completion
file_path
