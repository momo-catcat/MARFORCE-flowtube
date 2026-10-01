###### read the output profile file as an input for the next stage ######
import numpy as np


def read_output_profile(file_path,obj):
    with open(file_path, 'r') as f:
        raw_lines = [line.strip() for line in f if not line.startswith('#') and line.strip()]

    merged_arrays = []
    current_array = []

    for line in raw_lines:
        line = line.strip()
        if line.startswith('['):
            # Start a new array
            if current_array:
                # If there was a previous array being built, merge and add it
                merged_arrays.append(' '.join(current_array))
                current_array = []
            # Remove the opening bracket and add to current array
            current_array.append(line[1:])
        elif line.endswith(']'):
            # Remove the closing bracket and add to current array
            current_array.append(line[:-1])
            # Merge and add the completed array
            merged_arrays.append(' '.join(current_array))
            current_array = []
        else:
            # Middle line of an array
            current_array.append(line)

    # Handle case where the last array wasn't closed properly
    if current_array:
        merged_arrays.append(' '.join(current_array))

    # Now merged_arrays contains each complete array as a single string
    # for arr in merged_arrays:
    #     print(f"[{arr}]")
    import re
    data = []
    for arr_str in merged_arrays:
        # Fix the typo (replace '?' with '0')
        arr_str = arr_str.replace('?', '0')
        # Extract all numbers using regex
        numbers = re.findall(r'[-+]?\d*\.?\d+e[-+]?\d+', arr_str)
        # Convert to floats
        float_numbers = [float(num) for num in numbers]
        data.append(float_numbers)

    # Step 3: Create a pandas DataFrame
    # df = pd.DataFrame(data)
    data = np.array(data)
    data = data.reshape(obj.Rgrid,obj.Zgrid,obj.comp_num)

    return data