import os
import shutil

# Specify the directory where your files are located
directory = '/home/perera/Develop/ACESII/acespy/libso'

# List all files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.a'):
        # Build the new filename by replacing '.a' with '.so'
        new_filename = filename.replace('.a', '.so')

        # Construct the full file paths
        old_path = os.path.join(directory, filename)
        new_path = os.path.join(directory, new_filename)

        # Rename the file
        shutil.move(old_path, new_path)
        print(f'Renamed: {filename} to {new_filename}')

