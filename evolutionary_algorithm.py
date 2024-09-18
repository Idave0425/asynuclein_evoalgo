import random
import subprocess
import tempfile
import sys
import os

POSSIBLE_AA = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
POP_SIZE = 30
# Population Initialization
# NAC region has 35 AA, we choose ~10 before and after to make 55 AA molecule
# Call this function like initialize_population(30, 55)
def initialize_population(n, m):
  pop_i = [[random.choice(POSSIBLE_AA) for _ in range(m)] for _ in range(n)]
  return pop_i

# Fitness Function - binary addition to 15
#def calculate_specificity(individual):
#    # INSERT_YOUR_CODE
#
#    def run_mi_smart_filters(individual):
#        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as temp_fasta1, \
#             tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta'):
#            # Write the individual to both temporary fasta files
#            fasta_content = f">individual\n{''.join(individual)}\n"
#            temp_fasta1.write(fasta_content)
#            temp_fasta2.write(fasta_content)
#            temp_fasta1_path = temp_fasta1.name
#            temp_fasta2_path = temp_fasta2.name
#
#        # Define the command to run mi_smart_filters.py
#        command = [
#            sys.executable, 'pic_specificity_algorithm/mi_smart_filters.py',
#            '--fasta1', temp_fasta1_path,
#            '--fasta2', temp_fasta2_path,
#            '--dataname', 'output.csv',
#            '--imagename', 'heatmap.png',
#            '--protein1', 'protein1',
#            '--protein2', 'protein2',
#            '--modelSpecies', 'individual'
#        ]
#
#        # Run the command and capture the output
#        result = subprocess.run(command, capture_output=True, text=True)
#
#        # Clean up temporary files
#        os.remove(temp_fasta1_path)
#        os.remove(temp_fasta2_path)
#
#        if result.returncode != 0:
#            print(f"Error running mi_smart_filters.py: {result.stderr}")
#            return None
#
#        # Process the output CSV to extract the specificity score
#        specificity_score = 0
#        with open('output.csv', 'r') as csvfile:
#            for line in csvfile:
#                if line.startswith('individual'):
#                    specificity_score = float(line.split(',')[2])  # Assuming the score is in the 3rd column
#
#        return specificity_score

#    total = run_mi_smart_filters(individual)

#    return total

def calculate_specificity(individual):
    aa_ranking = {aa: rank for rank, aa in enumerate(sorted(set(individual)), start=1)}
    specificity_score = sum(aa_ranking[aa] for aa in individual)
    return specificity_score

# Function to remove x individuals with lowest scores
def highestscore_selection(population, scores, x):
    sorted_indices = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    updated_population = [population[i] for i in sorted_indices[:x]]
    return updated_population


# Function to generate new individuals near existing ones
def generate_new_individuals(population, pop_size, X):
    new_generation = population[:]
    while len(new_generation) < pop_size:
        parent = random.choice(population)
        new_individual = list(parent)
        indices_to_change = random.sample(range(len(parent)), X)
        for index in indices_to_change:
            new_individual[index] = random.choice(POSSIBLE_AA)  # Assuming POSSIBLE_AA is a list of amino acids
        new_generation.append(new_individual)  # Append the list of amino acids directly
    return new_generation

# Iterative Filtering and Modification
def optimize_molecule(num_gen, pop_size, aa_len, sel_str, mut_rate):
    pop_i = initialize_population(pop_size, aa_len)
    pop = pop_i[:]
    for i in range(num_gen):
        print(f"Generation {i}")
        print(f"population {pop}")
        pop_specificities = []
        for individual in pop:
            specificity_score = calculate_specificity(individual)
            pop_specificities.append(specificity_score)
        selected_indivuals = highestscore_selection(pop_i, pop_specificities, sel_str)
        new_generation = generate_new_individuals(selected_indivuals, POP_SIZE, mut_rate)
        pop = new_generation[:]
    return pop_specificities, pop



# Example usage:
target_individual = find_target_individual(10, 5)
print("Target Individual:", target_individual)