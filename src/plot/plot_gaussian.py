import numpy as np
import matplotlib.pyplot as plt
import re
import os  # Added for directory creation


def read_file(filename):
    """Read the content of the file and return it as a list of lines"""
    try:
        with open(filename, 'r') as file:
            content = file.readlines()
        return content
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None


def read_parameters(filename):
    """Read parameters from the parameter file"""
    try:
        with open(filename, 'r') as file:
            content = file.readlines()
            
        # First line contains the length
        length = float(content[0].strip())
        
        # Second line contains the positions
        positions = [float(x) for x in content[1].strip().split()]
        
        return length, positions
    except FileNotFoundError:
        print(f"Error: Parameter file '{filename}' not found.")
        return 10.0, [0.0]  # Default values if file not found
    except (IndexError, ValueError) as e:
        print(f"Error reading parameter file: {e}")
        return 10.0, [0.0]  # Default values if file format is incorrect


def parse_basis_set(content):
    """Parse the basis set information from the file content"""
    basis_info = {
        'element': None,
        's_functions': [],
        'p_functions': []
    }
    
    # Variables to track parsing state
    current_section = None
    num_primitives = 0
    num_contracted = 0
    
    i = 0
    while i < len(content):
        line = content[i].strip()
        
        # Skip empty lines
        if not line:
            i += 1
            continue
        
        # Check for element info
        if line.startswith('A '):
            basis_info['element'] = {'symbol': 'C', 'atomic_number': int(line.split()[1])}
            i += 1
            continue
        
        # Check for section headers
        if line.startswith('$ S-TYPE FUNCTIONS'):
            current_section = 's'
            i += 1
            # Next line contains counts
            counts = content[i].strip().split()
            num_primitives = int(counts[0])
            num_contracted = int(counts[1])
            
            # Initialize data structures for S functions
            exponents = []
            coefficients = [[] for _ in range(num_contracted)]
            
            # Parse the primitive data
            for j in range(num_primitives):
                i += 1
                values = re.findall(r'[-+]?\d*\.\d+|\d+', content[i])
                exponents.append(float(values[0]))
                for k in range(num_contracted):
                    if k+1 < len(values):
                        coefficients[k].append(float(values[k+1]))
            
            # Store the parsed functions
            for k in range(num_contracted):
                primitives = []
                for j in range(num_primitives):
                    if j < len(exponents) and k < len(coefficients) and j < len(coefficients[k]):
                        primitives.append({
                            'exponent': exponents[j],
                            'coefficient': coefficients[k][j]
                        })
                
                if primitives:
                    basis_info['s_functions'].append({
                        'type': 's',
                        'primitives': primitives
                    })
            
        elif line.startswith('$ P-TYPE FUNCTIONS'):
            current_section = 'p'
            i += 1
            # Next line contains counts
            counts = content[i].strip().split()
            num_primitives = int(counts[0])
            num_contracted = int(counts[1])
            
            # Initialize data structures for P functions
            exponents = []
            coefficients = [[] for _ in range(num_contracted)]
            
            # Parse the primitive data
            for j in range(num_primitives):
                i += 1
                values = re.findall(r'[-+]?\d*\.\d+|\d+', content[i])
                exponents.append(float(values[0]))
                for k in range(num_contracted):
                    if k+1 < len(values):
                        coefficients[k].append(float(values[k+1]))
            
            # Store the parsed functions
            for k in range(num_contracted):
                primitives = []
                for j in range(num_primitives):
                    if j < len(exponents) and k < len(coefficients) and j < len(coefficients[k]):
                        primitives.append({
                            'exponent': exponents[j],
                            'coefficient': coefficients[k][j]
                        })
                
                if primitives:
                    basis_info['p_functions'].append({
                        'type': 'p',
                        'primitives': primitives
                    })
        
        i += 1
    
    return basis_info


def compute_gaussian(x, exponent, coefficient, center=0.0):
    """Compute the value of a Gaussian primitive at x"""
    return coefficient * np.exp(-exponent * (x - center) ** 2)


def get_periodic_distance(x, center, L):
    """
    Compute the correct periodic distance between x and center in a box of length L
    Returns the proper signed distance considering periodicity
    
    This function is vectorized to handle array inputs properly
    """
    # Direct distance
    direct_dist = x - center
    
    # Initialize the result array with the direct distance
    result = direct_dist.copy() if isinstance(direct_dist, np.ndarray) else direct_dist
    
    # Calculate the periodically wrapped distance for each point
    
    # For scalar inputs
    if np.isscalar(direct_dist):
        if direct_dist < -L/2:  # If going right across boundary is shorter
            result = direct_dist + L
        elif direct_dist > L/2:  # If going left across boundary is shorter
            result = direct_dist - L
    # For array inputs
    else:
        # If going right across boundary is shorter
        mask_right = direct_dist < -L/2
        result[mask_right] = direct_dist[mask_right] + L
        
        # If going left across boundary is shorter
        mask_left = direct_dist > L/2
        result[mask_left] = direct_dist[mask_left] - L
    
    return result


def compute_periodic_gaussian(x, exponent, coefficient, center=0.0, L=10.0):
    """
    Compute the value of a periodic Gaussian primitive at x
    
    The periodic Gaussian has the form:
    g(x) = coefficient * exp(-exponent * d^2)
    
    where d is the minimum distance between x and center in a periodic box of length L
    """
    # Get the periodic distance (this is the signed distance considering periodicity)
    d_periodic = get_periodic_distance(x, center, L)
    
    # Compute the Gaussian using the squared distance (sign doesn't matter for exponential part)
    return coefficient * np.exp(-exponent * d_periodic**2)


def compute_contracted_gaussian(x, primitives, center=0.0, periodic=False, L=10.0):
    """Compute the value of a contracted Gaussian at x"""
    result = np.zeros_like(x, dtype=float)
    for primitive in primitives:
        if periodic:
            result += compute_periodic_gaussian(x, primitive['exponent'], primitive['coefficient'], center, L)
        else:
            result += compute_gaussian(x, primitive['exponent'], primitive['coefficient'], center)
    return result


def plot_gaussian_overlap(basis_info, length, positions, periodic=False, output_dir="plot_tmp"):
    """Plot and shade the overlap between Gaussian functions at different positions"""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a high-resolution grid for plotting
    x = np.linspace(0, length, 1000)
    
    # Need at least two positions to show overlap
    if len(positions) < 2:
        print("Need at least two positions to show overlap")
        return
    
    # Process s-type orbitals
    for i, s_func in enumerate(basis_info['s_functions']):
        plt.figure(figsize=(10, 8))
        
        # Create the orbital name
        orbital_name = f"{i+1}s"
        
        # Create a list to store all position values for later overlap calculation
        all_values = []
        
        # Calculate and plot values for each position
        for position in positions:
            values = compute_contracted_gaussian(
                x, s_func['primitives'], center=position, 
                periodic=periodic, L=length
            )
            plt.plot(x, values, '-', linewidth=2, label=f"Position {position}")
            all_values.append(values)
        
        # Shade overlap regions between consecutive pairs
        for j in range(len(positions)-1):
            values1 = all_values[j]
            values2 = all_values[j+1]
            
            # Find the minimum value at each point (the overlap)
            overlap = np.minimum(values1, values2)
            
            # Shade the overlap area without a label
            plt.fill_between(x, 0, overlap, color=f'C{j}', alpha=0.3)
        
        title_prefix = "Periodic" if periodic else "Regular"
        plt.title(f"{title_prefix} {orbital_name} Orbital with Shaded Overlap Regions", fontsize=16)
        plt.xlabel('Distance (Bohr)', fontsize=14)
        plt.ylabel('Amplitude', fontsize=14)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        
        # Save with a unique name and close the figure
        prefix = "periodic" if periodic else "regular"
        filename = os.path.join(output_dir, f'{prefix}_{orbital_name}_orbital.png')
        plt.savefig(filename, dpi=300)
        plt.close()
        print(f"Saved: {filename}")
    
    # Process p-type orbitals
    for i, p_func in enumerate(basis_info['p_functions']):
        plt.figure(figsize=(10, 8))
        
        # Create the orbital name
        orbital_name = f"{i+1}p"
        
        # Create a list to store all position values for later overlap calculation
        all_values = []
        
        # Calculate values for each position and normalize
        for position in positions:
            if periodic:
                # For each Gaussian primitive
                p_values = np.zeros_like(x, dtype=float)
                
                for primitive in p_func['primitives']:
                    exponent = primitive['exponent']
                    coefficient = primitive['coefficient']
                    
                    # Get the correct periodic distance (preserving sign)
                    d_periodic = get_periodic_distance(x, position, length)
                    
                    # Calculate the p-orbital: (x-xA) * exp(-α(x-xA)²)
                    gaussian_value = coefficient * np.exp(-exponent * d_periodic**2)
                    p_orbital_value = d_periodic * gaussian_value
                    
                    p_values += p_orbital_value
                
            else:
                # For regular p-type orbitals, compute the contracted gaussian
                values = compute_contracted_gaussian(
                    x, p_func['primitives'], center=position,
                    periodic=False
                )
                
                # Apply the p-orbital angular factor (x - center)
                p_values = values * (x - position)
            
            # Normalize to have a maximum absolute value of 1
            max_val = np.max(np.abs(p_values))
            if max_val > 0:
                p_values = p_values / max_val
            
            plt.plot(x, p_values, '-', linewidth=2, label=f"Position {position}")
            all_values.append(p_values)
        
        # Shade overlap regions between consecutive pairs
        for j in range(len(positions)-1):
            values1 = all_values[j]
            values2 = all_values[j+1]
            
            # For p-orbitals, we need to handle positive and negative values
            # Calculate overlap where both have the same sign
            same_sign_mask = np.sign(values1) == np.sign(values2)
            overlap = np.zeros_like(x)
            
            # Where they have the same sign, take the minimum absolute value and keep sign
            overlap = np.where(same_sign_mask, 
                              np.where(np.abs(values1) < np.abs(values2), values1, values2),
                              0)
            
            # Shade the positive overlap area without a label
            pos_mask = overlap > 0
            if np.any(pos_mask):
                plt.fill_between(x, 0, overlap, where=pos_mask, color=f'C{j}', alpha=0.3)
            
            # Shade the negative overlap area without a label
            neg_mask = overlap < 0
            if np.any(neg_mask):
                plt.fill_between(x, overlap, 0, where=neg_mask, color=f'C{j}', alpha=0.3)
        
        title_prefix = "Periodic" if periodic else "Regular"
        plt.title(f"{title_prefix} {orbital_name} Orbital with Shaded Overlap Regions (Normalized)", fontsize=16)
        plt.xlabel('Distance (Bohr)', fontsize=14)
        plt.ylabel('Normalized Amplitude', fontsize=14)
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        
        # Save with a unique name and close the figure
        prefix = "periodic" if periodic else "regular"
        filename = os.path.join(output_dir, f'{prefix}_{orbital_name}_orbital.png')
        plt.savefig(filename, dpi=300)
        plt.close()
        print(f"Saved: {filename}")

def plot_gaussian_pairs(basis_info, length, positions, periodic=False, output_dir="plot_tmp"):
    """Plot pairs of different Gaussian orbitals at different positions"""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a high-resolution grid for plotting
    x = np.linspace(0, length, 1000)
    
    # Need at least two positions to show pairs
    if len(positions) < 2:
        print("Need at least two positions to show pairs")
        return
    
    # Process s-type orbitals pairs (different orbital types)
    for i in range(len(basis_info['s_functions']) - 1):
        for j in range(i + 1, len(basis_info['s_functions'])):
            s_func1 = basis_info['s_functions'][i]
            s_func2 = basis_info['s_functions'][j]
            
            # For each pair of positions
            for pos_idx1 in range(len(positions) - 1):
                for pos_idx2 in range(pos_idx1 + 1, len(positions)):
                    pos1 = positions[pos_idx1]
                    pos2 = positions[pos_idx2]
                    
                    plt.figure(figsize=(10, 8))
                    
                    # Calculate values for orbital 1 at position 1
                    values1 = compute_contracted_gaussian(
                        x, s_func1['primitives'], center=pos1, 
                        periodic=periodic, L=length
                    )
                    
                    # Calculate values for orbital 2 at position 2
                    values2 = compute_contracted_gaussian(
                        x, s_func2['primitives'], center=pos2, 
                        periodic=periodic, L=length
                    )
                    
                    # Normalize both s-orbitals
                    max_val1 = np.max(np.abs(values1))
                    if max_val1 > 0:
                        values1 = values1 / max_val1
                    
                    max_val2 = np.max(np.abs(values2))
                    if max_val2 > 0:
                        values2 = values2 / max_val2
                    
                    # Plot both Gaussians
                    plt.plot(x, values1, 'b-', linewidth=2, label=f"{i+1}s at Position {pos1}")
                    plt.plot(x, values2, 'r-', linewidth=2, label=f"{j+1}s at Position {pos2}")
                    
                    # Find the minimum value at each point (the overlap)
                    overlap = np.minimum(values1, values2)
                    
                    # Shade the overlap area
                    plt.fill_between(x, 0, overlap, color='purple', alpha=0.3, label="Overlap")
                    
                    # Set title and labels
                    title_prefix = "Periodic" if periodic else "Regular"
                    plt.title(f"{title_prefix} Pair {i+1}s-{j+1}s Orbitals at Positions {pos1}-{pos2}", 
                             fontsize=16)
                    plt.xlabel('Distance (Bohr)', fontsize=14)
                    plt.ylabel('Normalized Amplitude', fontsize=14)  # Updated label
                    plt.grid(True)
                    plt.legend()
                    plt.tight_layout()
                    
                    # Save with a unique name and close the figure
                    prefix = "periodic" if periodic else "regular"
                    filename = os.path.join(output_dir, 
                                           f'{prefix}_pair_{i+1}s_{j+1}s_pos_{pos1}_{pos2}.png')
                    plt.savefig(filename, dpi=300)
                    plt.close()
                    print(f"Saved: {filename}")
    
    # Process p-type orbitals pairs (different orbital types)
    for i in range(len(basis_info['p_functions']) - 1):
        for j in range(i + 1, len(basis_info['p_functions'])):
            p_func1 = basis_info['p_functions'][i]
            p_func2 = basis_info['p_functions'][j]
            
            # For each pair of positions
            for pos_idx1 in range(len(positions) - 1):
                for pos_idx2 in range(pos_idx1 + 1, len(positions)):
                    pos1 = positions[pos_idx1]
                    pos2 = positions[pos_idx2]
                    
                    plt.figure(figsize=(10, 8))
                    
                    # Calculate p-orbital values for position 1
                    if periodic:
                        # For each Gaussian primitive
                        p_values1 = np.zeros_like(x, dtype=float)
                        
                        for primitive in p_func1['primitives']:
                            exponent = primitive['exponent']
                            coefficient = primitive['coefficient']
                            
                            # Get the correct periodic distance (preserving sign)
                            d_periodic = get_periodic_distance(x, pos1, length)
                            
                            # Calculate the p-orbital: (x-xA) * exp(-α(x-xA)²)
                            gaussian_value = coefficient * np.exp(-exponent * d_periodic**2)
                            p_orbital_value = d_periodic * gaussian_value
                            
                            p_values1 += p_orbital_value
                    else:
                        # For regular p-type orbitals
                        values = compute_contracted_gaussian(
                            x, p_func1['primitives'], center=pos1,
                            periodic=False
                        )
                        
                        # Apply the p-orbital angular factor (x - center)
                        p_values1 = values * (x - pos1)
                    
                    # Calculate p-orbital values for position 2
                    if periodic:
                        # For each Gaussian primitive
                        p_values2 = np.zeros_like(x, dtype=float)
                        
                        for primitive in p_func2['primitives']:
                            exponent = primitive['exponent']
                            coefficient = primitive['coefficient']
                            
                            # Get the correct periodic distance (preserving sign)
                            d_periodic = get_periodic_distance(x, pos2, length)
                            
                            # Calculate the p-orbital: (x-xA) * exp(-α(x-xA)²)
                            gaussian_value = coefficient * np.exp(-exponent * d_periodic**2)
                            p_orbital_value = d_periodic * gaussian_value
                            
                            p_values2 += p_orbital_value
                    else:
                        # For regular p-type orbitals
                        values = compute_contracted_gaussian(
                            x, p_func2['primitives'], center=pos2,
                            periodic=False
                        )
                        
                        # Apply the p-orbital angular factor (x - center)
                        p_values2 = values * (x - pos2)
                    
                    # Normalize both p-orbitals
                    max_val1 = np.max(np.abs(p_values1))
                    if max_val1 > 0:
                        p_values1 = p_values1 / max_val1
                    
                    max_val2 = np.max(np.abs(p_values2))
                    if max_val2 > 0:
                        p_values2 = p_values2 / max_val2
                    
                    # Plot both p-orbitals
                    plt.plot(x, p_values1, 'b-', linewidth=2, label=f"{i+1}p at Position {pos1}")
                    plt.plot(x, p_values2, 'r-', linewidth=2, label=f"{j+1}p at Position {pos2}")
                    
                    # Calculate overlap where both have the same sign
                    same_sign_mask = np.sign(p_values1) == np.sign(p_values2)
                    overlap = np.where(same_sign_mask, 
                                      np.where(np.abs(p_values1) < np.abs(p_values2), p_values1, p_values2),
                                      0)
                    
                    # Shade the positive overlap
                    pos_mask = overlap > 0
                    if np.any(pos_mask):
                        plt.fill_between(x, 0, overlap, where=pos_mask, color='purple', alpha=0.3)
                    
                    # Shade the negative overlap
                    neg_mask = overlap < 0
                    if np.any(neg_mask):
                        plt.fill_between(x, overlap, 0, where=neg_mask, color='purple', alpha=0.3)
                    
                    # Add a single overlap entry to the legend
                    plt.plot([], [], color='purple', alpha=0.3, label="Overlap")
                    
                    # Set title and labels
                    title_prefix = "Periodic" if periodic else "Regular"
                    plt.title(f"{title_prefix} Pair {i+1}p-{j+1}p Orbitals at Positions {pos1}-{pos2}", 
                             fontsize=16)
                    plt.xlabel('Distance (Bohr)', fontsize=14)
                    plt.ylabel('Normalized Amplitude', fontsize=14)
                    plt.grid(True)
                    plt.legend()
                    plt.tight_layout()
                    
                    # Save with a unique name and close the figure
                    prefix = "periodic" if periodic else "regular"
                    filename = os.path.join(output_dir, 
                                           f'{prefix}_pair_{i+1}p_{j+1}p_pos_{pos1}_{pos2}.png')
                    plt.savefig(filename, dpi=300)
                    plt.close()
                    print(f"Saved: {filename}")



def main():
    basis_filename = "./tmp/Basis_scratch"
    param_filename = "./tmp/parameter.dat"
    plot_dir = "plot_tmp"  # Directory for saving plots
    
    print(f"Reading basis file: {basis_filename}")
    print(f"Reading parameter file: {param_filename}")
    print(f"Plots will be saved to: {plot_dir}/")
    
    # Read the files
    content = read_file(basis_filename)
    if content is None:
        return
    
    # Read parameters
    length, positions = read_parameters(param_filename)
    print(f"Plot length: {length}")
    print(f"Gaussian positions: {positions}")
    
    # Parse the basis set information
    basis_info = parse_basis_set(content)
    
    # Print summary of basis functions
    print(f"\nBasis set for {basis_info['element']['symbol']} (Z={basis_info['element']['atomic_number']})")
    
    print("\nS-type Functions:")
    for i, func in enumerate(basis_info['s_functions']):
        print(f"  {i+1}s Orbital:")
        print(f"    Number of primitives: {len(func['primitives'])}")
        print("    Exponents and Coefficients:")
        for prim in func['primitives']:
            print(f"      {prim['exponent']:.10f}  {prim['coefficient']:.10f}")
    
    print("\nP-type Functions:")
    for i, func in enumerate(basis_info['p_functions']):
        print(f"  {i+1}p Orbital:")
        print(f"    Number of primitives: {len(func['primitives'])}")
        print("    Exponents and Coefficients:")
        for prim in func['primitives']:
            print(f"      {prim['exponent']:.10f}  {prim['coefficient']:.10f}")
    
    # Plot regular Gaussian overlap
    print("\nGenerating regular Gaussian overlap plots...")
    plot_gaussian_overlap(basis_info, length, positions, periodic=False, output_dir=plot_dir)
    
    # Plot periodic Gaussian overlap
    print("\nGenerating periodic Gaussian overlap plots...")
    plot_gaussian_overlap(basis_info, length, positions, periodic=True, output_dir=plot_dir)
    
    print(f"\nAll plots have been saved to {plot_dir}/ folder. No plots were displayed.")

    plot_gaussian_pairs(basis_info, length, positions, periodic=False, output_dir=plot_dir)

    plot_gaussian_pairs(basis_info, length, positions, periodic=True, output_dir=plot_dir)



if __name__ == "__main__":
    main()
