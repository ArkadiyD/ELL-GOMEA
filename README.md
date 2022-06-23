# ELL-GOMEA
Code repository for the article 
> "Arkadiy Dushatskiy, Tanja Alderliesten, and Peter A. N. Bosman. 2021. A Novel Approach to Designing Surrogate-assisted Genetic Algorithms by Combining Efficient Learning of Walsh Coefficients and Dependencies. ACM Trans. Evol. Learn. Optim. 1, 2, Article 5 (June 2021), 23 pages. https://doi.org/10.1145/3453141"

To run Efficient Linkage Learning (ELL) algorithm:
1. Compile GOMEA, ELL with using  ```makefile_gomea, makefile_ell``` files
2. Run GOMEA on the specific problem to get the baseline performance: ```python3 run_algorithms.py GOMEA [problem_index] [first_run] [last_run]``` 
3. Run **ELL** (the best performing algorithm) (or alternatively *PLL, LARSLL*): ```python3 run_algorithms.py [algorithm_name] [problem_index] [first_run] [last_run]```
4. Check out results in the corresponding folder with name ```results/[algorithm_name]/[problem_name]/[problem_size]/[run_index]``` The fitness of the elitist solution can be found in the ```elitistonly.dat``` file. The elitist solution itself can be found in the ```elitistsolution.dat```.

