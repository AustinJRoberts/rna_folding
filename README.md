# RNA Folding

In biology and chemistry,
the properties of a molecule are not solely determined by a set of atoms
but also by the shape of the molecule. 
In genetics, the shape of an RNA molecule is largely determined by how it bends back on itself. 
The sequence of A’s, U’s, G’s, and C’s that make up RNA has 
certain pairs are drawn together to form hydrogen bonds.
A sequence of several bonds in a row is called a stem,
and a stem provides sufficient force to keep the molecule folded together.
In nature, an RNA molecule will form some stems while avoiding others
in a manner that  minimizes the free energy of the system. 

This demo program takes am RNA sequence and applies a quadratic model in pursuit of the optimal stem configuration.

<p align = "center">

![Figure 1!](readme_imgs/Single_Stem.png "This is a title")
<p align = "center">
Fig.1 - An RNA sequence with a single stem of length 4.
</p>

Predicting the existence of stems is important to predicting the properties of the RNA molecule.
However, prediction is made more complicated by two important factors. 

First, stems are not allowed to overlap. 
A simple case of overlapping can be illustrated by Figure 1.
The stem expressed by the pink lines can be denoted with the tuple (2, 5, 14, 17),
where 2 is the index of the first G,
4 is the index of the U,
14 is the index of the A,
and 17 is the index of the last C.
This 4-tuple is mapped to a single variable.
However, the smaller stems (2, 4, 15, 17) and (3, 5, 15, 16) also need to be considered,
even though the optimal solution will not include them in this case.

Second, the intertwining phenomenon known as pseudoknots are less energetically favorable.
In Figure 2, we see an example of such a pseudoknot, 
where one side of a stem occurs in between the two sides of a different stem.
The use of a quadratic objective allows us to make pseudoknots less likely to occur in optimal solutions,
increasing overall accuracy.

<p align = "center">

![Figure 2](readme_imgs/pseudoknot2.png)
<p align = "center">
Fig.2 - A pseudoknot formed by a stem of length 3 and a stem of length 5.
</p>

This demo is loosely based on the work in [1],
which is in turn inspired by [2].

## Usage

To run the demo through a command line interface, type:

```bash
python RNA_folding.py
```

The demo prints the optimal stem configuration along with other relevant data.
It then saves a plot of the sequence and its bonds as `RNA_plot.png`.

### Optional parameters
Several optional parameters are accepted:

- `--path`: specifies path the input text file with RNA sequence information. 
- `--verbose`: if set to default value of 'True',
the program prints additional information about the model. 
- `--min-stem`: minimum length necessary for a stem to be considered.
- `--min-loop`: minimum number of neucleotides that must be present
in between the two sides of a stem for that stem to be considered. 
In the literature, this is termed a 'hairpin loop.'
- `-c`: used in th coefficient, *ck<sub>i</sub>k<sub>j</sub>*, 
applied to the quadratic pseudoknot terms.
Larger values make pseudoknots less likely.

As an example, to explicitly call the default values, type:
```bash
python RNA_folding.py --path RNA_text_files/TMGMV_UPD-PK1.txt --verbose True  --min-stem 3 --min-loop 2 -c 0.3 
```


## Problem Formulation

In predicting the stems of an RNA molecule, we build a quadratic model with three contributing factors. 

1. Variables correspond to potential stems, 
linearly weighted by the negative square of their length, *k*.

2. Potential pseudoknots correspond to quadratic terms 
with weight equal to the product of the two lengths 
times a positive parameter *c*.

3. Overlapping stems are disallowed, which is enforced by a constraint.

![objective](readme_imgs/objective.png) 

Subject to ![constraint](readme_imgs/constraint.png) if stems *i* and *j* overlap.

Here, each *x<sub>i</sub>* is a binary variable indicating the inclusion/exclusion of the *i<sup>th</sup>* stem.
Each constant *k<sub>i</sub>* is the length of said stem.
The indexing set *S* is the set of all pairs of stems that forma a pseudoknot.
Finally, *c* is a tunable parameter adjusting the impact of pseudonknots.
It is set to 0.3 by default.

This formulation is loosely based on [1].

In the printed solution, each stem is denoted by four numbers. 
The first two numbers correspond to the beginning and ending indices of the first side of the stem. 
Similarly, the last two numbers correspond to the beginning and ending indices of the second side of the stem.

## Code Overview

The implementation can be broken into three main parts
1. Preprocessing the RNA sequence to extract all possible stems, pseudoknots, and overlaps.
2. Building the model and sending it to a hybrid solver to find a solution.
3. Post-processing the solution 
to print appropriate information and create the plot.

A majority of the code is dedicated to step 1. 
Here, possible bonds are stored in a binary matrix,
and the matrix is searched for possible stems.
Possible stems (each corresponding to a decision variable) 
are stored in a dictionary structure to speed processing.

## Code Specifics

By default, the minimum stem length is set to 3. 
A stem of length 5 thus contains
two stems of length 4 and three stems of length 3 under inclusion.
The stem dictionary records the maximal stems (under inclusion) as keys,
where each key maps to a list of the associated stems weakly contained within the maximal stem.

No two stems contained in the same maximal key can both be in an optimal solution, 
so we treat them all as overlapping, regardless of if it is literally the case.
This particular case of overlapping is enforced through a one-hot constraint to improve solver performance.

We further use the stem dictionary structure 
to avoid comparing all combinations of stems when searching for pseudoknots and overlaps.

Plotting uses a randomized process to find a reasonable layout. 
For this reason, the plot will change in successive runs, 
even if the solution does not. 

## References

[1] Fox DM, MacDermaid CM, Schreij AM, Zwierzyna M, Walker RC. 
RNA folding using quantum computers. 
bioRxiv; 2021. DOI: 10.1101/2021.05.27.446060.

[2] Kai, Zhang, et al. 
"An efficient simulated annealing algorithm for the RNA secondary structure prediction with Pseudoknots." 
BMC genomics 20.13 (2019): 1-13.