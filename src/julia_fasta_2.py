import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import logging
from collections import Counter
from multiprocessing import Pool, cpu_count
from functools import partial

def create_directories(directories):
    """
    Creates the necessary directories if they don't exist.

    Parameters:
        directories (list): A list of directory paths to create.
    """
    for directory in directories:
        try:
            os.makedirs(directory, exist_ok=True)
            print(f"Directory '{directory}' is ready.")
        except Exception as e:
            print(f"Error creating directory '{directory}': {e}")

def configure_logging(log_dir, log_file):
    """
    Configures the log settings.

    Parameters:
        log_dir (str): Directory where the log file will be stored.
        log_file (str): Name of the log file.
    """
    log_path = os.path.join(log_dir, log_file)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info("Logging configured successfully.")

def log_message(message):
    """
    Message to both console and log file.

    Parameters:
        message (str): The message to log.
    """
    print(message)
    logging.info(message)

def read_fasta(file_path):
    """
    Reads a FASTA file and extracts sequences.

    Parameters:
        file_path (str): Path to the FASTA file.

    Returns:
        list: A list of tuples (header, sequence).
    """
    sequences = []
    try:
        with open(file_path, 'r') as f:
            header = None
            seq = ''
            for line in f:
                if line.startswith('>'):
                    if header and seq:
                        sequences.append((header, seq))
                        seq = ''
                    header = line[1:].strip()
                else:
                    seq += line.strip()
            if header and seq:
                sequences.append((header, seq))
    except FileNotFoundError:
        log_message(f"Error: File {file_path} not found.")
    except Exception as e:
        log_message(f"Error reading FASTA file: {e}")
    return sequences

def calculate_gc_content(sequence):
    """
    Calculates the GC content of a sequence.

    Parameters:
        sequence (str): The nucleotide or amino acid sequence.

    Returns:
        float: GC content (%).
    """
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    total_bases = len([base for base in sequence.upper() if base in ['A', 'T', 'C', 'G']])
    if total_bases == 0:
        return 0.0
    return (gc_count / total_bases) * 100

def generate_k_mers(sequence, k):
    """
    Generates non-overlapping k-mers from a given sequence.

    Parameters:
        sequence (str): The nucleotide or amino acid sequence.
        k (int): Length of each k-mer.

    Returns:
        list: List of k-mers.
    """
    k_mers = [sequence[i:i+k] for i in range(0, len(sequence) - k + 1, k)]
    return k_mers

def k_mer_to_complex(k_mer):
    """
    Maps a single k-mer to a complex number.

    Parameters:
        k_mer (str): A k-length substring of the sequence.

    Returns:
        complex: A complex number derived from the k-mer.
    """
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    total_real, total_imag, valid_bases = 0, 0, 0
    for char in k_mer.upper():
        if char in purines:
            total_real += 1
            total_imag += 1
        elif char in pyrimidines:
            total_real -= 1
            total_imag += 1
        else:
            continue  # Skip non-standard bases
        valid_bases += 1
    if valid_bases == 0:
        return complex(0, 0)
    avg_real = total_real / valid_bases
    avg_imag = total_imag / valid_bases

    scale_factor = 0.5  # Adjusted
    real_part = avg_real * scale_factor
    imag_part = avg_imag * scale_factor

    return complex(real_part, imag_part)

def validate_and_adjust_c(c, min_mag=0.4, max_mag=1.0):
    """
    Validates and adjusts the complex constant c to lie within [min_mag, max_mag].

    Parameters:
        c (complex): The complex constant to validate.
        min_mag (float): Minimum magnitude.
        max_mag (float): Maximum magnitude.

    Returns:
        complex: Adjusted complex constant.
    """
    magnitude = abs(c)
    if magnitude == 0:
        return complex(0.4, 0.0)
    angle = np.angle(c)
    magnitude = np.clip(magnitude, min_mag, max_mag)
    return magnitude * complex(np.cos(angle), np.sin(angle))

def seq_to_complex(sequence, k, gc_content, seq_length):
    """
    Maps a sequence to a complex number using k-mer mapping, GC content and sequence length.

    Parameters:
        sequence (str): The nucleotide or amino acid sequence.
        k (int): Length of each k-mer.
        gc_content (float): GC content percentage of the sequence.
        seq_length (int): Length of the sequence.

    Returns:
        complex: Complex number derived from the sequence.
    """
    k_mers = generate_k_mers(sequence, k)
    if not k_mers:
        return complex(0, 0)

    kmer_counts = Counter(k_mers)
    total_kmers = sum(kmer_counts.values())

    weighted_real, weighted_imag = 0, 0
    for k_mer, count in kmer_counts.items():
        c = k_mer_to_complex(k_mer)
        weight = count / total_kmers  # Frequency-based weighting
        weighted_real += c.real * weight
        weighted_imag += c.imag * weight

    avg_c = complex(weighted_real, weighted_imag)

    # Integrate GC content and sequence length into c
    gc_normalized = gc_content / 100  # Scale GC content to [0,1]
    length_normalized = min(seq_length / 1000, 1.0)  # Scale sequence length to [0,1], cap at 1000

    # Adjust c based on GC content and sequence length
    adjusted_real = avg_c.real * (1 + gc_normalized)
    adjusted_imag = avg_c.imag * (1 + length_normalized)

    adjusted_c = complex(adjusted_real, adjusted_imag)

    # Validate and adjust c to ensure it lies within the desired magnitude range
    c = validate_and_adjust_c(adjusted_c)

    return c

def generate_julia(c, width, height, iterations):
    """
    Generates a Julia set image for a given complex number c.

    Parameters:
        c (complex): The complex parameter for the Julia set.
        width (int): Image width in pixels.
        height (int): Image height in pixels.
        iterations (int): Number of iterations for the algorithm.

    Returns:
        numpy.ndarray: An array representing the Julia set image.
    """
    x = np.linspace(-1.5, 1.5, width)
    y = np.linspace(-1.5, 1.5, height)
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y

    div_time = np.zeros(Z.shape, dtype=int)

    for _ in range(iterations):
        mask = np.abs(Z) <= 2
        Z[mask] = Z[mask] ** 2 + c
        div_time += mask

    return div_time

def normalize_julia(julia):
    """
    Normalizes the Julia set array for better visualization.

    Parameters:
        julia (numpy.ndarray): The Julia set array.

    Returns:
        numpy.ndarray: Normalized Julia set array.
    """
    julia_normalized = np.log(julia + 1)  # Enhance contrast
    julia_normalized = (julia_normalized - julia_normalized.min()) / (julia_normalized.max() - julia_normalized.min())
    return julia_normalized

def plot_julia(julia_normalized, output_file, zoom, pan_x, pan_y, cmap='inferno'):
    """
    Plots and saves the Julia set image.

    Parameters:
        julia_normalized (numpy.ndarray): Normalized Julia set array.
        output_file (str): Path to save the image.
        zoom (float): Zoom level.
        pan_x (float): Pan offset on the x-axis.
        pan_y (float): Pan offset on the y-axis.
        cmap (str): Colormap for the image.
    """
    plt.figure(figsize=(12, 12), dpi=300)
    plt.axis('off')
    zoom_scaled = 1.5 / zoom
    plt.xlim(-zoom_scaled + pan_x, zoom_scaled + pan_x)
    plt.ylim(-zoom_scaled + pan_y, zoom_scaled + pan_y)
    plt.imshow(
        julia_normalized,
        cmap=cmap,
        extent=(
            -zoom_scaled + pan_x,
            zoom_scaled + pan_x,
            -zoom_scaled + pan_y,
            zoom_scaled + pan_y
        )
    )
    try:
        plt.savefig(output_file, bbox_inches='tight', pad_inches=0)
    except Exception as e:
        log_message(f"Error saving image {output_file}: {e}")
    finally:
        plt.close()

def process_sequence(idx, header, sequence, args):
    """
    Process a single sequence to generate and save its Julia set.

    Parameters:
        idx (int): Sequence index.
        header (str): Sequence header.
        sequence (str): Nucleotide or amino acid sequence.
        args: Parsed command-line arguments.
    """
    try:
        gc_content = calculate_gc_content(sequence)
        seq_length = len(sequence)
        c = seq_to_complex(sequence, args.kmer, gc_content, seq_length)
        log_message(f"[Sequence {idx}] Header: {header}")
        log_message(f"[Sequence {idx}] Length: {seq_length}")
        log_message(f"[Sequence {idx}] GC Content: {gc_content:.2f}%")
        log_message(f"[Sequence {idx}] Complex Constant c = {c}")
        julia = generate_julia(c, args.width, args.height, args.iterations)
        julia_normalized = normalize_julia(julia)
        output_file = os.path.join(args.output, f'julia_{idx}.png')
        plot_julia(julia_normalized, output_file, args.zoom, args.pan_x, args.pan_y, args.cmap)
        log_message(f"[Sequence {idx}] Saved Julia set image to {output_file}")
    except Exception as e:
        log_message(f"[Sequence {idx}] Error processing sequence: {e}")

def main():
    """
    Main function.
    """
    directories = ["logs", "plots", "results", "output/images"]
    create_directories(directories)

    configure_logging("logs", "julia_set_generator.log")

    # Parameters
    parser = argparse.ArgumentParser(
        description='Generate Julia sets from a FASTA file with k-mer mapping, GC content, and sequence length integration.'
    )
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('-o', '--output', help='Output directory', default='output/images')
    parser.add_argument('-w', '--width', type=int, help='Image width in pixels', default=800)
    parser.add_argument('-e', '--height', type=int, help='Image height in pixels', default=800)
    parser.add_argument('-i', '--iterations', type=int, help='Number of iterations', default=300)
    parser.add_argument('--zoom', type=float, help='Zoom level (default: 1.0)', default=1.0)
    parser.add_argument('--pan_x', type=float, help='Pan offset in x-axis (default: 0.0)', default=0.0)
    parser.add_argument('--pan_y', type=float, help='Pan offset in y-axis (default: 0.0)', default=0.0)
    parser.add_argument('-k', '--kmer', type=int, help='k-mer length for mapping (default: 4)', default=4)
    parser.add_argument('--cmap', type=str, help='Colormap for Julia sets (default: inferno)', default='inferno')
    parser.add_argument('--processes', type=int, help='Number of parallel processes (default: number of CPU cores)', default=cpu_count())
    args = parser.parse_args()

    # Logging
    log_message(f"Reading sequences from {args.input}")
    sequences = read_fasta(args.input)
    if not sequences:
        log_message("No sequences found. Exiting.")
        return

    if not os.path.exists(args.output):
        try:
            os.makedirs(args.output)
            log_message(f"Created output directory at {args.output}")
        except Exception as e:
            log_message(f"Error creating output directory: {e}")
            return

    log_message(f"Starting processing with {args.processes} parallel processes.")

    pool = Pool(processes=args.processes)
    process_func = partial(process_sequence, args=args)
    tasks = [(idx + 1, header, seq) for idx, (header, seq) in enumerate(sequences)]
    try:
        pool.starmap(process_func, tasks)
    except KeyboardInterrupt:
        log_message("Processing interrupted by user.")
    except Exception as e:
        log_message(f"An error occurred during multiprocessing: {e}")
    finally:
        pool.close()
        pool.join()
        log_message("Processing completed.")

if __name__ == '__main__':
    main()
