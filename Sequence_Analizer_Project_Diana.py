
"""
Author: Diana Castano
Date:  07th Jun 2024.
This program is designed for analyzing DNA and protein sequences. 
The Graphic user interface is used to execute the analysis on provided DNA and protein sequence files, respectively. 
Select the files to be analised and click on analize sequences. 
"""


import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def read_sequences_from_file(file_path):
    """
    Reads sequences from a given file, removing unwanted lines/identifiers.
    
    Args:
    file_path (str): The path to the input file containing sequences.

    Returns:
    str: A string of concatenated sequences with identifiers removed.
    """    
    sequences = []
    with open(file_path, 'r') as fin:
        sequence = ''
        for line in fin:
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                sequence = ''
            else:
                sequence += line.strip()
        if sequence:
            sequences.append(sequence)
    return ''.join(sequences).upper()

def count_nucleotides(seq):
    """
    Counts the occurrences of each nucleotide (A, T, G, C, N) in a sequence. N means non identified 
    
    Args:
    seq (str): The DNA sequence to count nucleotides in.

    Returns:
    dict: A dictionary with the count of each nucleotide.
    """    
    count_nuc = {"A":0, "T":0, "G":0, "C":0, "N":0}
    for nucleotide in seq:
        if nucleotide in count_nuc:
            count_nuc[nucleotide] += 1
    return count_nuc

def plot_nucleotide_counts(counts):
    """
        Plots the count of each nucleotide in a DNA sequence.
        
        Args:
        counts (dict): A dictionary containing the count of each nucleotide.
    
    """       
    nucleotides = list(counts.keys())
    values = list(counts.values())
    axes = plt.axes()
    axes.bar(nucleotides, values, color=['blue', 'green', 'red', 'purple', 'orange'])
    axes.set_title('Nucleotide Counts')
    axes.set_xlabel('Nucleotide')
    axes.set_ylabel('Count')
    plt.show()

def read_protein_sequences(file_path):
    """
    Reads protein sequences from a FASTA file.
    
    Args:
    file_path (str): The path to the input FASTA file containing protein sequences.

    Returns:
    list: A list of protein sequences.
    """    
    protein_sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        protein_sequences.append(str(record.seq))
    return protein_sequences

def calculate_sequence_lengths(sequences):
    """
    Calculates the lengths of protein sequences.
    
    Args:
    sequences (list): A list of protein sequences.

    Returns:
    list: A list of lengths corresponding to each protein sequence.
    """    
    return [len(seq) for seq in sequences]

def analyze_amino_acid_composition(sequences):
    """
    Analyzes the amino acid composition of protein sequences.
    
    Args:
    sequences (list): A list of protein sequences.

    Returns:
    list: A list of dictionaries with the count of each amino acid for each protein sequence.
    """    
    return [ProteinAnalysis(seq).count_amino_acids() for seq in sequences]

def plot_sequence_lengths(lengths):
    """
    Plots the distribution of sequence lengths.
    
    Args:
    lengths (list): A list of sequence lengths.
    """    
    plt.hist(lengths, bins=20, color='skyblue', edgecolor='black')
    plt.title('Distribution of Sequence Lengths')
    plt.xlabel('Sequence Length')
    plt.ylabel('Frequency')
    plt.show()

def analyze_protein_sequence(protein_sequence):
    """
    Analyze the given protein sequence.
    
    Args:
    protein_sequence (str): The protein sequence to analyze.
    
    Returns:
    dict: A dictionary containing the count of each amino acid.
    """    
    protein_analysis = ProteinAnalysis(protein_sequence)
    return protein_analysis.count_amino_acids()

def plot_amino_acid_composition(amino_acid_composition):
    """
    Plot the amino acid composition using colors.
    
    Args:
    amino_acid_composition (dict): A dictionary containing the count of each amino acid.
    """    
    amino_acids = list(amino_acid_composition.keys())
    counts = list(amino_acid_composition.values())
    colors = plt.cm.tab20(np.linspace(0, 1, len(amino_acids)))
    plt.figure(figsize=(12, 6))
    plt.bar(amino_acids, counts, color=colors)
    plt.title('Amino Acid Composition')
    plt.xlabel('Amino Acid')
    plt.ylabel('Count')
    plt.show()


# GUI Part
class SequenceAnalyzerApp(tk.Tk):
    """
       A GUI application for analyzing DNA and protein sequences.
   
       Attributes:
           title (str): The title of the application window.
           geometry (str): The initial dimensions of the application window.
           file_paths (list): A list to store the paths of selected sequence files.
           protein_sequence (str): The protein sequence read from the file.
   
       Methods:
           __init__: Initializes the SequenceAnalyzerApp class and creates UI elements.
           select_files: Opens a file dialog to select sequence files for analysis.
           analyze_sequences: Analyzes the selected sequence files and displays results.
           read_protein_sequence_from_file: Reads a protein sequence from a file.
    """
    def __init__(self):
        """
        Initializes the SequenceAnalyzerApp class and creates UI elements.
        """        
        super().__init__()
        self.title("Sequence Analyzer")
        self.geometry("800x600")
        
        self.file_paths = []
        self.protein_sequence = ""
        
        # Create UI elements
        self.label = tk.Label(self, text="Select sequence files for analysis")
        self.label.pack(pady=10)
        
        self.select_files_btn = tk.Button(self, text="Select Files", command=self.select_files)
        self.select_files_btn.pack(pady=10)
        
        self.analyze_btn = tk.Button(self, text="Analyze Sequences", command=self.analyze_sequences)
        self.analyze_btn.pack(pady=10)
        
        self.result_text = scrolledtext.ScrolledText(self, width=100, height=20, wrap=tk.WORD)
        self.result_text.pack(pady=10)
        
    def select_files(self):
        """
        Opens a file dialog to select sequence files for analysis.
        """        
        self.file_paths = filedialog.askopenfilenames(filetypes=[("FASTA files", "*.fna *.fasta"), ("All files", "*.*")])
        if self.file_paths:
            messagebox.showinfo("Files Selected", f"{len(self.file_paths)} files selected")
        
    def analyze_sequences(self):
        """
        Analyzes the selected sequence files and displays results.
        """        
        if not self.file_paths:
            messagebox.showwarning("No Files", "No files selected for analysis")
            return
        
        self.result_text.delete(1.0, tk.END)
        
        for file_path in self.file_paths:
            self.result_text.insert(tk.END, f"Results for file: {file_path}\n")
            
            if file_path.endswith(".fna"):
                sequence = read_sequences_from_file(file_path)
                result = count_nucleotides(sequence)
                self.result_text.insert(tk.END, f"Nucleotide Counts:\n{result}\n")
                self.result_text.insert(tk.END, f"Total nucleotides counted: {sum(result.values())}\n\n")
                plot_nucleotide_counts(result)
                
            elif file_path.endswith(".fasta"):
                protein_sequences = read_protein_sequences(file_path)
                lengths = calculate_sequence_lengths(protein_sequences)
                plot_sequence_lengths(lengths)
                
                for i, sequence in enumerate(protein_sequences):
                    self.result_text.insert(tk.END, f"Sequence {i+1} - Amino Acid Composition:\n")
                    amino_acid_composition = analyze_protein_sequence(sequence)
                    self.result_text.insert(tk.END, f"{amino_acid_composition}\n\n")
                    plot_amino_acid_composition(amino_acid_composition)
            else:
                messagebox.showwarning("Invalid File", f"Unsupported file type: {file_path}")
                
    def read_protein_sequence_from_file(self, file_path):
        """
        Reads a protein sequence from a file.
        
        Args:
        file_path (str): The path to the input file containing the protein sequence.
        """
        try:
            with open(file_path, 'r') as file:
                self.protein_sequence = file.read().strip()
                messagebox.showinfo("Sequence Loaded", "Protein sequence successfully loaded.")
        except Exception as e:
            messagebox.showerror("Error", f"Error reading file: {str(e)}")


# Run the GUI application
app = SequenceAnalyzerApp()
app.mainloop()
