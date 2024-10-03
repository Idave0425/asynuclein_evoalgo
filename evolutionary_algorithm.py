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

def calculate_specificity(individual):
    """
    Calculate the specificity score of an individual amino acid sequence.

    The specificity score is determined by assigning a rank to each unique amino acid
    in the sequence based on its alphabetical order. The score is the sum of these ranks
    for all amino acids in the sequence. Higher scores indicate sequences with more 
    diverse amino acids.

    Args:
        individual (list): A list of amino acids representing an individual sequence.

    Returns:
        int: The specificity score of the individual sequence.
    """
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
        selected_individuals = highestscore_selection(pop, pop_specificities, sel_str)
        new_generation = generate_new_individuals(selected_individuals, pop_size, mut_rate)
        pop = new_generation[:]
    return pop_specificities, pop



# Example usage:
target_individual = find_target_individual(10, 5)
print("Target Individual:", target_individual)