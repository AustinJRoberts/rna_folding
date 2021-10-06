import os
from RNA_folding import *
import neal
import dimod
from dwave.system import LeapHybridSampler


def find_overlaps(stem_dict):
    """ Uses a stem dictionary to find all overlapping stems that are not contained in a common maximal stem.

    Args:
        stem_dict: Dictionary with maximal stems as keys and list of weakly contained substems as values.

    Returns: List of pairs of overlapping stems.
    """

    overlaps = []
    for stem1, stem2 in combinations(stem_dict.keys(), 2):
        # Check maximal stems first.
        if check_overlap(stem1, stem2):
            # If maximal stems overlap, compare list of smaller stems.
            overlaps.extend([(substem1, substem2) for substem1, substem2 in product(stem_dict[stem1], stem_dict[stem2])
                             if check_overlap(substem1, substem2)])

    return overlaps


def build_bqm(stem_dict, lam=None):
    """ Builds a binary quadratic model from a stem dictionary.

    Args:
        stem_dict: Dictionary with maximal stems as keys and list of weakly contained substems as values.
        lam: Float coefficient used as 'Legrangian' penatly coefficient for overlaps.
        c: Float coefficient used to penalize the existance of pseudoknots in solutions.

    Returns: Dimod BinaryQuadraticModel
    """

    if lam is None:
        lam = max(stem[1]-stem[0] + 1 for stem in stem_dict.keys())**2 + 3

    # Create linear coefficients the prioritize inclusion of long stems.
    linear_coeffs = {stem: -1*(stem[1] - stem[0] + 1)**2 for sublist in stem_dict.values() for stem in sublist}

    # Create constraints for overlapping and and substem containment.
    quadratic_coeffs = {stems: lam for stems in find_overlaps(stem_dict)}

    for substems in stem_dict.values():
        quadratic_coeffs.update({(stem1, stem2): lam for (stem1, stem2) in (combinations(substems, 2))})

    # Add penalty to discourage pseudoknots.
    quadratic_coeffs.update(pseudoknot_terms(stem_dict))

    # Build BQM from dictionaries.
    bqm = dimod.BinaryQuadraticModel(linear_coeffs, quadratic_coeffs, 'BINARY')

    return bqm


def process_bqm_solution(samples, verbose=True):
    """ Processes samples from solution.

    Prints information about the best feasible solution and returns a list of stems contained in solution.

    Args:
        samples: dimod sample element.
        verbose: Boolean indicating if function should print additional information.

    Returns: List of stems included in optimal solution, encoded as 4-tuples.
    """
    sample_set = samples.aggregate().data()

    feasible = False
    for sample in sample_set:
        bonded_stems = [stem for stem, val in sample[0].items() if val == 1]

        # Check feasibility of sample.
        if [stems for stems in combinations(bonded_stems, 2) if check_overlap(*stems)]:
            print('Warning: Ignoring best solution with infeasible values.')
            pass
        else:
            best_energy = sample[1]
            feasible = True
            break

    if not feasible:
        print('\nWarning! All solutions infeasible due to overlapping stems. Consider changing parameters.')
        return None

    print('\n# stems in best solution:', len(bonded_stems))
    print('Stems in best solution:', *bonded_stems)

    if verbose:
        print('\n# variables (stems):', len(sample[0].keys()))

        # find psuedoknots using product instead of combinations allows for short asymmetric checks.
        pseudoknots = [(stem1, stem2) for [stem1, stem2] in product(bonded_stems, bonded_stems)
                       if stem1[1] < stem2[0] and stem2[1] < stem1[2] and stem1[3] < stem2[2]]

        print('\n# pseudoknots in best solution:', len(pseudoknots))
        if pseudoknots:
            print('Pseudoknots:', *pseudoknots)

        print('\nBest Energy:', best_energy)

    return bonded_stems


if __name__ == "__main__":
    file = 'RNA_text_files/TMGMV_UPD-PK1.txt'  # min_stem 3, 6
    # file = 'RNA_text_files/NGF-L6.txt' # min_stem 4, 19 vs 18?
    # file = 'RNA_text_files/STMV_UPD2-PK1.txt' # min_stem 3, 59 stems
    # file = 'RNA_text_files/NC_008516.txt'  # min_stem 3, 177
    # file = 'RNA_text_files/hiv.txt'
    file = 'RNA_text_files/simple.txt'
    file = 'RNA_text_files/simple_pseudo.txt'  # min_stem 4, 19 vs 18?

    this_folder = os.path.dirname(os.path.abspath(__file__))
    path_to_file = os.path.join(this_folder, file)

    print('\nPulling data from:', path_to_file)

    stem_dict = make_stem_dict(text_to_matrix(path_to_file))

    bqm = build_bqm(stem_dict)
    sampler = neal.SimulatedAnnealingSampler()
    sample_set = sampler.sample(bqm, num_reads=20)

    stems = process_bqm_solution(sample_set)
    make_plot(path_to_file, stems)