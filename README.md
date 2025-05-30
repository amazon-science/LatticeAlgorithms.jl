## LatticeAlgorithms.jl

This package contains lattice algorithms that were used in the paper [Closest lattice point decoding for multimode Gottesman-Kitaev-Preskill codes](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.4.040334). The package contains several [lattice reduction algorithms](https://www.ant.uni-bremen.de/sixcms/media.php/102/10740/SPM_2011_Wuebben.pdf), such as [Lenstra-Lenstra-Lovász](https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm) and [Korkine-Zolotarev](https://en.wikipedia.org/wiki/Korkine%E2%80%93Zolotarev_lattice_basis_reduction_algorithm) algorithms, and a [search algorithm](https://publications.lib.chalmers.se/records/fulltext/14990/local_14990.pdf) for solving the [closest point problem](https://en.wikipedia.org/wiki/Lattice_problem#Closest_vector_problem_(CVP)) and the [shortest vector problem](https://en.wikipedia.org/wiki/Lattice_problem#Shortest_vector_problem_(SVP)). For the Gottesman-Kitaev-Preskill (GKP) codes, the package includes  the $D_n$ lattice and two types of repetition-GKP codes, which can be decoded efficiently from a lattice perspective. The data and code for the paper can be found in the folder `examples/papers/Closest_lattice_point_decoding_for_multimode_GKP_codes`.

This package also contains several algorithms that were used in the paper [Exploring the quantum capacity of a Gaussian random displacement channel using
Gottesman-Kitaev-Preskill codes and maximum likelihood decoding](https://arxiv.org/abs/2411.04277), including an efficient and exact maximum likelihood decoder (MLD) for surface-square GKP codes, and a tensor-network decoder to approximate the MLD for generic concatenated-GKP codes. The latter is built on top of the decoder in [SweepContractor.jl](https://github.com/chubbc/SweepContractor.jl). The data and code for the paper can be found in the folder `examples/papers/Exploring_the_quantum_capacity_of_a_Gaussian_random_displacement_channel_using_GKP_codes_and_maximum_likelihood_decoding`. 

This package also contains several algorithms that were used in the paper *Approximate maximum likelihood decoding with $K$ minimum weight matchings*  (to appear soon). In this paper, we introduce a novel algorithm to approximate the optimal maximum likelihood decoding strategey via finding $K$ minimum weight matchings from the decoding graph. The data and plots for the paper can be found in the folder `examples/papers/Approximate_maximum_likelihood_decoding_with_K_minimum_weight_matchings`. 



## Security

See [CONTRIBUTING](CONTRIBUTING.md#security-issue-notifications) for more information.

## Citation
If you find our paper or codebase useful, please consider citing us as:
```
@article{PRXQuantum.4.040334,
  title = {Closest Lattice Point Decoding for Multimode Gottesman-Kitaev-Preskill Codes},
  author = {Lin, Mao and Chamberland, Christopher and Noh, Kyungjoo},
  journal = {PRX Quantum},
  volume = {4},
  issue = {4},
  pages = {040334},
  numpages = {36},
  year = {2023},
  month = {Dec},
  publisher = {American Physical Society},
  doi = {10.1103/PRXQuantum.4.040334},
  url = {https://link.aps.org/doi/10.1103/PRXQuantum.4.040334}
}
```

## Examples

We provide some examples for using the package. The location to the source code of a function can be look up in the `src/LatticeAlgorithms.jl` file. For example, the function "closest_point" is exported before the line "include("lattice_algorithms.jl")", indicating that the function is in the file `src/lattice_algorithms.jl`. More tutorials can be found in the `examples/tutorials` folder. 

**Example 1**: Finding the closest point for a random lattice
```
using LatticeAlgorithms
n = 10 # Dimension of the lattice
M = rand(n, n) # A lattice generated by a random matrix
x = rand(n) # A random input vector
y = closest_point(x, M) # The closest lattice point to x
```
Finding the closest point lies in the core of solving other lattice problems, including shortest lattice vector problem, and finding the relevant vectors for a given lattice. The demonstrations of solving these problems can be found in the folder `examples/tutorials`. 

**Example 2**: Finding the closest point for root lattices
```
using LatticeAlgorithms
n = 20 # Dimension of the lattice

M = Dn(n)
x = rand(n)

@time y1 = closest_point(x, M) # 0.001310 seconds (18.21 k allocations: 633.391 KiB)
@time y2 = closest_point_Dn(x) # 0.000014 seconds (3 allocations: 672 bytes)
@assert y1 ≈ y2 # true
```
Finding the closest point for [root lattices](http://neilsloane.com/doc/Me83.pdf) can be done efficiently, in contrast to general lattices. In the above example, we demonstrate this fact with the $D_n$ lattice. 


**Example 3**: Lattice reduction
```
using LatticeAlgorithms
n = 10 # Dimension of the lattice
M = rand(n, n) # A lattice generated by a random matrix
lll_basis = lll(M)
@assert lll_basis.T * lll_basis.L * lll_basis.Q ≈ M # true
```
In the above example, we perform the LLL reduction to the given matrix. The output ```lll_basis``` contains three matrices such that ```T*L*Q=M```. Here ```T``` is a unimodular matrix, ```Q``` is an orthogonal matrix and ```L``` is the LLL basis that satisfies the LLL criteria. A given matrix can be checked if it satisfies the LLL criteria via
```
islllreduced(lll_basis.L) # true
```
Similarly the given matrix can be KZ reduced via
```
kz_basis = kz(M)
@assert kz_basis.T * kz_basis.L * kz_basis.Q ≈ M # true
```

**Example 4**: Finding the distances of Gottesman-Kitaev-Preskill (GKP) codes
```
using LatticeAlgorithms
M = surface_code_M(3)
distance(M)
```
In the above example, we find the code distances of a surface-square GKP code, which is defined as the Euclidean length of the shortest operators. To find the distances of X, Y and Z logical operators, we can use ```distances(M)```. More utilities regarding GKP codes, including canonization of lattice generator, finding logical operators and others can be found in the file ```src/gkp.jl```. 


**Example 5**: Exact maximum likelihood decoding (MLD) for the unrotated surface code on a square lattice
```
using LatticeAlgorithms
d = 25
n = d^2 + (d-1)^2

ϵs = 5/100 * ones(n)
prob_I, prob_X = bsv_unrotated_surface_code_qubit(ϵs, zeros(Int, n))

@assert (prob_I - log10(1.78283e-27)) < 1e-9 # true
@assert (prob_X - log10(5.58438e-57)) < 1e-9 # true
```
In the above example, we find the coset probability for a d=25 unrotated surface code with error rate 5%. The results are consistent with those given in [this paper](https://arxiv.org/pdf/1405.4883.pdf) by Bravyi, Suchara and Vargo (BSV, see TABEL I in page 16). The function `bsv_unrotated_surface_code_qubit` relies on a core subroutine `bsv` which is demonstrated below. More details for decoding rotated the (rotated) surface-square GKP code using the exact MLD can be found [here](https://github.com/amazon-science/LatticeAlgorithms.jl/blob/main/examples/papers/Exploring_the_quantum_capacity_of_a_Gaussian_random_displacement_channel_using_GKP_codes_and_maximum_likelihood_decoding/Get_data_bsv_surf_square.ipynb).


**Example 6**: Brute force approach decoding and further demonstration of BSV's exact MLD
```
using LatticeAlgorithms
d = 3
n = d^2 + (d-1)^2

stab_generators = unrotated_surface_code_Z_stabilizers(d)
stabilizers = get_stabilizer_group_from_generators(collect(values(stab_generators)))
ws = 0.6 * randn(n)
ws = abs.(ws) # The weight needs to be positive values

# Brute force approach to calculate the coset probability
expected_coset_prob = sum([prod(ws[stab]) for stab in stabilizers])

# BSV's exact method to calculate the log10 of the coset probability
coset_prob_bsv = bsv(ws)
@assert 10^value ≈ expected_value # true
```
In the above example, we find the coset probability for a d=3 unrotated surface code with two approaches, the brute force approach using the definition and BSV's approach. More details for decoding [[5-1-3]]-hexagonal GKP code using the brute force approach can be found [here](https://github.com/amazon-science/LatticeAlgorithms.jl/blob/main/examples/papers/Exploring_the_quantum_capacity_of_a_Gaussian_random_displacement_channel_using_GKP_codes_and_maximum_likelihood_decoding/get_data_mld_513_hex_code.ipynb).

**Example 7**: Tensor network decoding to approximate MLD for color-hexagonal GKP code
```
using LatticeAlgorithms

# Define the relevant quantities for a d=3 color-hexagonal GKP code
d = 3
n = triangular_color_code_num_qubits(d)
stab_generators = triangular_color_code_stabilizers(d)
full_stabilizers = get_stabilizer_group_from_generators(collect(values(stab_generators)))
S = [2 1; 0 sqrt(3)] / (12)^(1/4)
X = √π * S * [1, 0]
Z = √π * S * [0, 1]
bigX = vcat([X for _ in 1 : num_qubits]...) # The X logical operator
bigZ = vcat([Z for _ in 1 : num_qubits]...) # The X logical operator

σ = 0.6
ws = abs.(σ * randn(2n))

# Tensor network
χ = 64 # bond dimension

# Get the template for carrying out tensor network decoding
# The template can be repeatedly used for different syndrome
TN, indices = tn_template_color_hex(d)
lstar2, prob_I2, prob_X2, prob_Y2, prob_Z2, _, _, _, _ = tn_color_hex(ws, σ, TN, indices, bigZ, bigX, χ)
```
In the above example, we find the coset probability for a d=3 color-hexagonal GKP code with a variant of tensor-network decoder, which is built on top of the decoder in [SweepContractor.jl](https://github.com/chubbc/SweepContractor.jl). More details for decoding color-hexagonal GKP code using the tensor-network approach can be found [here](https://github.com/amazon-science/LatticeAlgorithms.jl/blob/main/examples/papers/Exploring_the_quantum_capacity_of_a_Gaussian_random_displacement_channel_using_GKP_codes_and_maximum_likelihood_decoding/get_data_tn_color_hex.ipynb).



## License

This project is licensed under the Apache-2.0 License.
