# Copyright 2021 D-Wave Systems
# Based on the paper 'RNA folding using quantum computers’
# Fox DM, MacDermaid CM, Schreij AM, Zwierzyna M, Walker RC.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
from itertools import product, combinations
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import dimod
from dwave.system import LeapHybridCQMSampler

"""
This file contains functions called  to help in the building and processing of models related to RNA folding.
The models are built and solved in a different file.
"""


def text_to_matrix(file_name, min_loop=2):
    """ Reads properly formatted RNA text file and returns a matrix of possible hydrogen bonding pairs.

    Args:
        file_name: String path to text file.
        min_loop: Integer minimum number of nucleotides between two bonding nucleotides.

    Returns: Numpy matrix of 0's and 1's, where 1 represents a possible bonding pair.
    """

    # Requires text file of rna data written in same format as examples.
    with open(file_name) as f:
        rna = "".join(("".join(line.split()[1:]) for line in f.readlines())).lower()

    # Create a dictionary of all indices where each nucleotide occurs.
    # Mote: t vs u is mainly preference of database. While u is in rna, t is in dna.
    index_dict = {'a': [], 'c': [], 'g': [], 't': [], 'u': []}

    # create a dictionary giving list of indices for each nucleotide.
    for i, nucleotide in enumerate(rna.lower()):
        index_dict[nucleotide].append(i)

    # List of possible hydrogen bonds for stems.
    # Recall that 't' is just a stand in for 'u'.
    hydrogen_bonds = [('a', 't'), ('a', 'u'), ('c', 'g'), ('g', 't'), ('g', 'u')]

    # Create a 0/1 matrix indicated where bonds may occur.
    bond_matrix = np.zeros((len(rna), len(rna)), dtype=bool)
    for pair in hydrogen_bonds:
        for bond in product(index_dict[pair[0]], index_dict[pair[1]]):
            if abs(bond[0] - bond[1]) > min_loop:
                bond_matrix[min(bond), max(bond)] = 1

    return bond_matrix


def make_stem_dict(bond_matrix, min_stem=3, min_loop=2):
    """ Takes a matrix of potential hydrogen binding pairs and returns a dictionary of possible stems.

    Each key is a maximal stem (under inclusion) whose value pair is a list of stems weakly contained in the key.
    Recording stems in this manner allows for faster computations.

    Args:
        bond_matrix: Numpy matrix of 0's and 1's, where 1 represents a possible bonding pair.
        min_stem: Integer minimum number of nucleotides between two bonding nucleotides.
        min_loop: Integer minimum number of nucleotides between two bonding nucleotides.

    Returns: Dictionary of all possible stems with maximal stems as keys.
    """

    stem_dict = {}
    n = bond_matrix.shape[0]

    # Iterate through matrix looking for possible stems.
    for i in range(n - (2 * min_stem + min_loop)):
        for j in range(i + 2 * min_stem + min_loop - 1, n):
            if bond_matrix[i, j]:
                k = 1
                while bond_matrix[i + k, j - k]:
                    bond_matrix[i + k, j - k] = False
                    k += 1

                if k >= min_stem:
                    # A 4-tuple is used to represent the stem.
                    stem_dict[(i, i + k - 1, j - k + 1, j)] = []

    # Iterate through all substems weakly contained in a maximal stem under inclusion.
    for stem in stem_dict.keys():
        stem_dict[stem].extend([(stem[0] + i, stem[0] + k, stem[3] - k, stem[3] - i)
                                for i in range(stem[1] - stem[0] - min_stem + 2)
                                for k in range(i + min_stem - 1, stem[1] - stem[0] + 1)])

    return stem_dict


def check_overlap(stem1, stem2):
    """ Checks if 2 stems use any of the same nucleotides.

    Args:
        stem1: 4-tuple containing stem information.
        stem2: 4-tuple containing stem information.

    Returns: Boolean indicating if the two stems overlap.
    """

    if type(stem1) == str or type(stem2) == str:
        return False

    # Check if any endpoints of stem2 overlap with stem1.
    for val in stem2:
        if stem1[0] <= val <= stem1[1] or stem1[2] <= val <= stem1[3]:
            return True
    # Check if endpoints of stem1 overlap with stem2.
    # Do not need to check all stem1 endpoints.
    for val in stem1[1:3]:
        if stem2[0] <= val <= stem2[1] or stem2[2] <= val <= stem2[3]:
            return True

    return False


def pseudoknot_terms(stem_dict, min_stem=3, c=0.3):
    """ Creates a dictionary with all possible pseudoknots as keys and appropriate penalty as as value pair.

    The penalty is the parameter c times the product of the lengths of the two stems in the knot.

    Args:
        stem_dict: Dictionary with maximal stems as keys and list of weakly contained substems as values.
        min_stem: Integer smallest number of consecutive bonds to be considered a stem.
        c: Float parameter factor of the penalty on pseudoknots.

    Returns: Dictionary with all possible pseudoknots as keys and appropriate penalty as as value pair.
    """

    pseudos = {}
    # Look within all pairs of maximal stems for possible psuedoknots.
    for stem1, stem2 in product(stem_dict.keys(), stem_dict.keys()):
        # Using product instead of combinations allows for short asymmetric checks.
        if stem1[0] + 2 * min_stem < stem2[1] and stem1[2] + 2 * min_stem < stem2[3]:
            pseudos.update({(substem1, substem2): c * (substem1[1] - substem1[0]) * (substem2[1] - substem2[0])
                            for substem1, substem2
                            in product(stem_dict[stem1], stem_dict[stem2])
                            if substem1[1] < substem2[0] and substem2[1] < substem1[2] and substem1[3] < substem2[2]})
    return pseudos


def make_plot(file, stems, fig_name='RNA_plot'):
    """ Produces graph plot and saves as .png file.

    Args:
        file: String text file name containing RNA information.
        stems: List of stems in solution, encoded as 4-tuples.
        fig_name: String name of file created to save figure. ".png" is added automatically

    Returns:None
    """

    # Read rna file for length and labels.
    with open(file) as f:
        rna = "".join(("".join(line.split()[1:]) for line in f.readlines())).lower()

    # Create graph with edges from RNA sequence and stems. Nodes are temporarily labeled by integers.
    G = nx.Graph()
    rna_edges = [(i, i + 1) for i in range(len(rna) - 1)]
    stem_edges = [(stem[0] + i, stem[3] - i) for stem in stems for i in range(stem[1] - stem[0] + 1)]
    G.add_edges_from(rna_edges + stem_edges)

    # Assign each nucleotide to a color.
    color_map = []
    for node in rna:
        if node == 'g':
            color_map.append('tab:red')
        elif node == 'c':
            color_map.append('tab:green')
        elif node == 'a':
            color_map.append('y')
        else:
            color_map.append('tab:blue')

    options = {"edgecolors": "tab:gray", "node_size": 250, "alpha": 0.8}
    pos = nx.spring_layout(G, iterations=5000)  # max(3000, 125 * len(rna)))
    nx.draw_networkx_nodes(G, pos, node_color=color_map, **options)

    labels = {i: rna[i].upper() for i in range(len(rna))}
    nx.draw_networkx_labels(G, pos, labels, font_size=10, font_color="whitesmoke")

    nx.draw_networkx_edges(G, pos, edgelist=rna_edges, width=3.0, alpha=0.5)
    nx.draw_networkx_edges(G, pos, edgelist=stem_edges, width=4.5, alpha=0.8, edge_color='tab:pink')

    plt.savefig(fig_name + '.png')


def build_cqm(stem_dict, min_stem=3):
    """ Creates a Constrained Binary Model to optimize most likely stems from a dictionary of possible stems.

    Args:
        stem_dict: Dictionary with maximal stems as keys and list of weakly contained substems as values.
        min_stem: Integer minimum number of nucleotides between two bonding nucleotides.

    Returns: Constrained Binary Model.
    """

    # Create linear coefficients the prioritize inclusion of long stems.
    linear_coeffs = {stem: -1 * (stem[1] - stem[0] + 1) ** 2 for sublist in stem_dict.values() for stem in sublist}

    # Create constraints for overlapping and and substem containment.
    quadratic_coeffs = pseudoknot_terms(stem_dict)

    bqm = dimod.BinaryQuadraticModel(linear_coeffs, quadratic_coeffs, 'BINARY')

    cqm = dimod.ConstrainedQuadraticModel()
    cqm.set_objective(bqm)

    # Add constraint disallowing overlapping substsems included in same maximal stem.
    for stem, substems in stem_dict.items():
        if len(substems) > 1:
            # Add the variable for all zeros case in one-hot constraint
            zeros = 'Null:' + str(stem)
            cqm.add_variable(zeros, 'BINARY')
            cqm.add_discrete(substems + [zeros], stem)

    overlaps = []
    for stem1, stem2 in combinations(stem_dict.keys(), 2):
        # Check maximal stems first.
        if check_overlap(stem1, stem2):
            # If maximal stems overlap, compare list of smaller stems.
            for stem_pair in product(stem_dict[stem1], stem_dict[stem2]):
                if check_overlap(stem_pair[0], stem_pair[1]):
                    cqm.add_constraint(dimod.quicksum([dimod.Binary(stem) for stem in stem_pair]) <= 1)

    # Add constraint disallowing other overlapping stems.
    # for stem_pairs in find_overlaps(stem_dict):
    #     cqm.add_constraint(dimod.quicksum([dimod.Binary(stem) for stem in stem_pairs]) <= 1)

    return cqm


def process_cqm_solution(sample_set, verbose=True):
    """ Processes samples from solution and prints relevant information.

    Prints information about the best feasible solution and returns a list of stems contained in solution.
    Returns solution as a list of stems rather than a binary string.

    Args:
        sample_set: dimod sampleset element.
        verbose: Boolean indicating if function should print additional information.

    Returns: List of stems included in optimal solution, encoded as 4-tuples.
    """

    infeasible = True
    for sample in sample_set.data():
        if sample.is_feasible:
            solution = sample
            infeasible = False
            break

    # Check for feasible solution.
    if infeasible:
        print('\nWarning! All solutions infeasible. You may need to try again.')
        return None

    print('Best Energy:', solution[1])

    # extract stems with a positive indicator variable.
    bonded_stems = [stem for stem, val in solution[0].items() if val == 1 and type(stem) == tuple]

    print('\n# stems in best solution:', len(bonded_stems))
    print('Stems in best solution:', *bonded_stems)

    if verbose:
        print('\n# variables (stems):', len(solution[0].keys()))

        # find psuedoknots using product instead of combinations allows for short asymmetric checks.
        pseudoknots = [(stem1, stem2) for [stem1, stem2] in product(bonded_stems, bonded_stems)
                       if stem1[1] < stem2[0] and stem2[1] < stem1[2] and stem1[3] < stem2[2]]

        print('\n# pseudoknots in best solution:', len(pseudoknots))
        if pseudoknots:
            print('Pseudoknots:', *pseudoknots)

    return bonded_stems


if __name__ == "__main__":
    file = 'RNA_text_files/TMGMV_UPD-PK1.txt'
    # file = 'RNA_text_files/NGF-L6.txt'
    # file = 'RNA_text_files/STMV_UPD2-PK1.txt'
    # file = 'RNA_text_files/NC_008516.txt'
    # file = 'RNA_text_files/hiv.txt'
    # file = 'RNA_text_files/simple.txt'

    this_folder = os.path.dirname(os.path.abspath(__file__))
    path_to_file = os.path.join(this_folder, file)

    print('\nPreprocessing data from:', path_to_file)

    stem_dict = make_stem_dict(text_to_matrix(path_to_file))
    cqm = build_cqm(stem_dict)

    print('Connecting to Solver...')
    sampler = LeapHybridCQMSampler()
    print('Finding Solution...')
    sample_set = sampler.sample_cqm(cqm)

    print('Processing solution...')
    stems = process_cqm_solution(sample_set)

    make_plot(path_to_file, stems)