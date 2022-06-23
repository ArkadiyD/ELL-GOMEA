#pragma once

#include <cmath>
#include <algorithm>
#include <numeric>
using namespace std;

#include "utils.hpp"
#include "problems.hpp"

struct Solution
{
	vector<bool> genotype;
	double fitness;
	
	Solution(){};
	Solution(int number_of_variables)
	{
		genotype.resize(number_of_variables);
	}

	Solution& operator= (const Solution &solution)
	{
	    genotype.resize(solution.genotype.size());
	    copy(solution.genotype.begin(), solution.genotype.end(), genotype.begin());
	    fitness = solution.fitness;
	    return *this;
	}
};


class GBOOptimizer
{
	int number_of_variables, algorithm;
	long long time_limit;
	Problem *problem_instance;
	vector<vector<int> > walsh_coefficients_by_index;
	vector<vector<int> > walsh_coefficients;
	vector<double> walsh_coefficients_values;
	unordered_set <vector<int>, hash_int_vector > linkage_pairs;
	unordered_map<vector<int>, double, hash_int_vector>  wCoef;
	double VTR;
	string folder;
	Solution elitist;
	long long elitist_number_of_evaluations;
	long long elitist_time;	
	bool elitist_solution_written_before = 0;
	long long timestamp_start;
	double number_of_evaluations = 0;
	
	void write_elitist_to_file();
	void check_elitist_and_vtr(Solution &solution);
	double installedProblemEvaluation(Solution &solution);
	double installedProblemEvaluationPartial(Solution &solution, vector<int> &touched_genes_indices, Solution &previous_solution);
	void dfs(int start_vertex, vector<vector<int> > &graph, vector<int> &labels, int component);
	void find_connected_components(vector<vector<int> > &graph, vector<vector<int> > &components);
	void partition_crossover(Solution &solution1, Solution &solution2, Solution &child);
	void hamming_ball_hill_climber(Solution &solution);
	void perturb(Solution &solution, double alpha);

	void optimize_surrogate_with_drils(double alpha);
	void optimize_surrogate_with_hirels();
	//void optimize_surrogate_with_population_px();

public:
	GBOOptimizer(int number_of_variables_, int algorithm_, char *folder_, unordered_map<vector<int>, double, hash_int_vector> &wCoef, unordered_set <vector<int>, hash_int_vector > &linkage_pairs);
	
	void optimize(long long time_limit_, double VTR_);
};

