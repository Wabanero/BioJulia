import os
import random
from Bio import SeqRecord, Seq
import numpy as np
import matplotlib.pyplot as plt
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import mean_squared_error as mse
from collections import defaultdict

# Create necessary directories
directories = ["logs", "plots", "results"]
for directory in directories:
    if not os.path.exists(directory):
        os.makedirs(directory)

# Function to log messages and print to console
def log_message(message):
    print(message)
    with open("logs/log.txt", "a") as log_file:
        log_file.write(message + "\n")

# Function to generate random 256 bp genomic sequences
def generate_random_sequence(length=256):
    return ''.join(random.choices('ATCG', k=length))

# Function to convert nucleotide to complex number
def nucleotide_to_complex(nucleotide):
    mapping = {'A': 1 + 1j, 'T': 1 - 1j, 'C': -1 + 1j, 'G': -1 - 1j}
    return mapping.get(nucleotide, 0 + 0j)

# Function to convert gene to complex sequence
def gene_to_complex_sequence(gene):
    return [nucleotide_to_complex(nuc) for nuc in gene]

# Function to normalize complex sequence
def normalize_complex_sequence(complex_sequence):
    max_value = max(abs(c) for c in complex_sequence)
    if max_value != 0:
        return [c / max_value for c in complex_sequence]
    return complex_sequence

# Function to adjust complex constant
def adjust_complex_constant(c):
    # Ensuring the complex constant is within a range that generates detailed fractals
    if abs(c) < 0.2:
        c *= 5
    elif abs(c) > 2:
        c /= 5
    return c

# Function to iterate Julia set
def iterate_julia(c, z, max_iter):
    for n in range(max_iter):
        if abs(z) > 2:
            return n
        z = z*z + c
    return max_iter

# Function to generate Julia set
def generate_julia_set(gene_sequence, width, height, zoom, move_x, move_y, max_iter):
    julia = np.zeros((width, height))
    gene_sequence = normalize_complex_sequence(gene_sequence)
    c = sum(gene_sequence) / len(gene_sequence)
    c = adjust_complex_constant(c)
    log_message(f"Adjusted complex constant c for Julia set: {c}")
    for x in range(width):
        for y in range(height):
            zx = 1.5 * (x - width / 2) / (0.5 * zoom * width) + move_x
            zy = 1.0 * (y - height / 2) / (0.5 * zoom * height) + move_y
            z = zx + zy * 1j
            julia[x, y] = iterate_julia(c, z, max_iter)
    return julia

# Function to compare two Julia sets
def compare_julia_sets(julia_set1, julia_set2):
    ssim_index = ssim(julia_set1, julia_set2)
    mse_value = mse(julia_set1, julia_set2)
    return ssim_index, mse_value

# Function to save unique nucleotide differences as FASTA files
def save_differences_as_fasta(differences_dict, sequences_dict):
    for (species1, idx1, species2, idx2), differences in differences_dict.items():
        gene_id1 = sequences_dict[species1][idx1].id
        gene_id2 = sequences_dict[species2][idx2].id
        filename = f"results/{species1.replace(' ', '_')}_vs_{species2.replace(' ', '_')}_differences.fasta"
        with open(filename, "w") as output_handle:
            for pos, nuc1, nuc2 in differences:
                header = f">{gene_id1}_pos{pos}_vs_{gene_id2}_pos{pos}"
                sequence = f"{nuc1}->{nuc2}"
                output_handle.write(f"{header}\n{sequence}\n")
        log_message(f"Saved differences between {species1} and {species2} to {filename}")

# Generate two random sequences
random_sequences = {
    "Random_Seq_1": [SeqRecord.SeqRecord(Seq.Seq(generate_random_sequence()), id="Random_Seq_1")],
    "Random_Seq_2": [SeqRecord.SeqRecord(Seq.Seq(generate_random_sequence()), id="Random_Seq_2")]
}

# Print the generated sequences
log_message(f"Generated Sequence 1 (Random_Seq_1): {random_sequences['Random_Seq_1'][0].seq}")
log_message(f"Generated Sequence 2 (Random_Seq_2): {random_sequences['Random_Seq_2'][0].seq}")

# Dictionary to store Julia sets for comparison
julia_sets = {}
sequences_dict = random_sequences

# Generate Julia sets for the random sequences
for species in random_sequences:
    julia_sets[species] = []
    for idx, record in enumerate(random_sequences[species]):
        plot_file = f"plots/{species.replace(' ', '_')}_sequence_{idx+1}.png"
        gene_sequence = record.seq
        complex_sequence = gene_to_complex_sequence(gene_sequence)
        if not complex_sequence:
            continue
        julia_set = generate_julia_set(complex_sequence, 1000, 1000, 1, 0, 0, 1000)
        julia_sets[species].append(julia_set)

        # Plot the Julia set with inverted color scheme
        plt.figure(figsize=(10, 10), dpi=300)
        plt.imshow(julia_set, extent=[-1.5, 1.5, -1, 1], cmap='inferno', origin='lower')
        plt.gca().invert_yaxis()
        plt.title(f"Julia Set for {species} (Sequence {idx+1})")
        plt.xlabel("Re(z)")
        plt.ylabel("Im(z)")
        plt.colorbar()
        plt.savefig(plot_file, dpi=300)
        plt.close()
        log_message(f"Saved Julia set plot for {species} sequence {idx+1} to {plot_file}.")

# Compare Julia sets between the random sequences
comparison_results = defaultdict(list)
for species1 in random_sequences:
    for species2 in random_sequences:
        if species1 != species2:
            for idx1, julia1 in enumerate(julia_sets[species1]):
                for idx2, julia2 in enumerate(julia_sets[species2]):
                    ssim_index, mse_value = compare_julia_sets(julia1, julia2)
                    comparison_results[(species1, species2)].append((idx1, idx2, ssim_index, mse_value))
                    log_message(f"Comparison of {species1} Sequence {idx1+1} with {species2} Sequence {idx2+1}:")
                    log_message(f"SSIM: {ssim_index:.4f}, MSE: {mse_value:.4f}")

# Identify nucleotide variations responsible for unique Julia sets
unique_differences = defaultdict(list)
for (species1, species2), comparisons in comparison_results.items():
    for idx1, idx2, ssim_index, mse_value in comparisons:
        if ssim_index < 0.9:  # Arbitrary threshold for significant difference
            seq1 = sequences_dict[species1][idx1]
            seq2 = sequences_dict[species2][idx2]
            differences = find_sequence_differences(seq1, seq2)
            unique_differences[(species1, idx1, species2, idx2)] = differences

# Save unique nucleotide differences as FASTA files
save_differences_as_fasta(unique_differences, sequences_dict)

log_message("Analysis complete. All files saved.")
