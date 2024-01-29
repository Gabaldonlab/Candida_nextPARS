import os
import matplotlib.pyplot as plt

def parse_dotbracket(dotbracket):
    # Extract pairs from dot-bracket notation
    pairs = []
    stack = []
    for i, char in enumerate(dotbracket):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                pairs.append((stack.pop(), i))
    return pairs  # Reverse the list to have lower index first

def calculate_long_range_interactions(pairs, length, percent_threshold=25):
    # Calculate threshold based on percentage of sequence length
    threshold = length * percent_threshold / 100
    # Count long-range interactions
    long_range_count = sum(1 for pair in pairs if pair[1] - pair[0] > threshold)
    

    return long_range_count / len(pairs) * 100

def plot_long_range_interactions(file_path):
    with open(file_path, 'r') as file:
        sequences = file.read().split('\n>')

    filenames = []
    percentages = []
    for sequence in sequences:
        if not sequence:
            continue  # Skip empty sequences
        lines = sequence.split('\n')
        filename = lines[0].split()[-1]
        filenames.append(filename)
        dotbracket = lines[2]
        pairs = parse_dotbracket(dotbracket)
        sequence_length = len(lines[1])
        percentage = calculate_long_range_interactions(pairs, sequence_length)
        percentages.append(percentage)

	# Get the PDF file name from file_path
    pdf_filename = "output_2023/" + os.path.splitext(os.path.basename(file_path))[0] + '.pdf'
    
    # Plotting
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(filenames, percentages, marker='o', linestyle='None')
    ax.set_xlabel('File Name')
    ax.set_ylabel('Percentage of Long-Range Interactions')
    ax.set_xticklabels(filenames, rotation=45, ha='right')
    ax.set_ylim(0, 100)  # Set y-axis limit from 0 to 100
    plt.subplots_adjust(bottom=0.2)
    plt.savefig(pdf_filename, dpi=300)
    
	# Get the PDF file name from file_path
    pdf_filename = "output_2023/" + os.path.splitext(os.path.basename(file_path))[0] + '.line.pdf'
    
    # Plotting
    fig, ax = plt.subplots(figsize=(5, 6))
    ax.plot(filenames, percentages, marker='o')
    ax.set_xlabel('File Name')
    ax.set_ylabel('Percentage of Long-Range Interactions')
    ax.set_xticklabels(filenames, rotation=45, ha='right')
    ax.set_ylim(0, 100)  # Set y-axis limit from 0 to 100
    plt.subplots_adjust(bottom=0.2)
    plt.savefig(pdf_filename, dpi=300)
    # ~ plt.show()

# Replace 'your_file_path_here' with the actual path to your file
file_path = 'your_file_path_here'

plot_long_range_interactions(file_path)



