import re
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


### IMPORT FILES

opts = ["2", "5", "8"]
dirs = ["inflam+minrad", "inflam+minrad+radcap", "inflam+rad_range"]
objs = ["Floor", "Ceiling", "Range"]

k = 0
legend_elements = []
fig, ax = plt.subplots()

for opt, dir, obj in zip(opts, dirs, objs):
    
    bestJhist = os.path.join('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/fig4', dir , opt, 'results', 'BestJhist.dat')
    log = os.path.join('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/fig4', dir, opt, 'log.txt')
    Jhist = os.path.join('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/fig4', dir, opt, 'results', 'Jhist.dat')
    Hhist = os.path.join('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/fig4', dir, opt, 'results', 'Hhist.dat')

    ### EXTRACT FROM LOG

    # Define the regular expression pattern to match the desired numbers
    pattern = r"Number of function evaluations till now :   (\d+)"

    # Initialize an empty list to store the extracted numbers

    numbers = []

    # Open and read the file
    with open(log, 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                # Extract the number and append it to the list
                numbers.append(int(match.group(1)))

    # Print the array of extracted numbers
    print(numbers)

    ### EXTRACT FROM JHIST

    # Initialize an empty list to store the extracted numbers
    jhist_numbers = []

    # Open and read the file
    with open(bestJhist, 'r') as file:
        for line in file:
            # Split the line by whitespace and convert the first element to float
            number = float(line.split()[0])
            jhist_numbers.append(number)

    ### FINDING CONSTRAINT VALUES

    Hhist_numbers = []

    # Loop through the elements in 'nums'
    for num in jhist_numbers:
        # Find the line number of 'num' in Jhist.dat
        with open(Jhist, 'r') as file:
            line_number = 1  # Initialize line number to 1
            for line in file:
                if float(line.strip()) == num:  # Assuming the values in Jhist.dat are float numbers
                    break
                line_number += 1

        # Get the corresponding value from Hhist.dat
        with open(Hhist, 'r') as file:
            Hhist_val = None
            for i, line in enumerate(file):
                if i + 1 == line_number:  # Match line numbers
                    Hhist_val = float(line.strip())
                    break

        if Hhist_val is not None:
            Hhist_numbers.append(Hhist_val)

    # 'values_from_Hhist' now contains the values from Hhist.dat corresponding to the elements in 'nums'
    print(Hhist_numbers)

    ### PLOTTING

    # Assuming you already have the 'numbers' array from the previous code

    # Plotting
    cmap = plt.get_cmap('hsv')
    colors = cmap(np.linspace(0,0.67,3))
    legend_elements.append(Line2D([0], [0], marker='o', color=colors[k], label=obj))

    plt.plot(numbers, jhist_numbers, color='k', linestyle='--')
    plt.scatter(numbers, jhist_numbers, color='white')
    for x, y, c in zip(numbers, jhist_numbers, Hhist_numbers):
        if c < 0:
            plt.plot(x, y, marker='o', color=colors[k], label=obj)
        else:
            plt.plot(x, y, marker='o', color=colors[k], fillstyle='none')
    
    k += 1

ax.legend(handles=legend_elements, loc='upper right')
plt.xlabel('Number of Function Calls')
plt.ylabel('Objective Function Value')
plt.title("Figure 4")
plt.savefig("fig4_objs.pdf")    