#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>
#include <functional>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <string>
using namespace std;

#include "utils.hpp"
#include "problems.hpp"
#include "gbo.hpp"

double EPS = 1e-6;
double R2_THRESHOLD;
int USE_HECKENDORN_ALGORITHM, SKIP_ALREADY_EVALUATED, OPTIMIZATION_ALGORITHM;
long long RANDOM_SEED;
int number_of_evaluations = 0,
  number_of_variables,
  problem_index,
  MAX_ORDER,
  max_evaluated_solutions = 10000000;
long        timestamp_start,                                                    /* The time stamp in milliseconds for when the program was started. */
            timestamp_finish;                                         /* The time stamp in milliseconds for when the algorithm was started (after problem initialization). */
char folder[1000];
long long gomea_time = 0;
double best_r2 = -1e+308;
int total_number_of_probes = 0;

unordered_map< vector<bool>, double, hash_vector > evaluated_solutions;
vector<double> evaluated_archive;
map<vector<int>, double> all_zeros_probes;

int holdout_population_size = number_of_variables;
long int timeLimitMilliseconds, timeLimitMillisecondsOptimization;

vector<vector<bool> > population;
vector<double> fitness_values;
vector<vector<bool> > holdout_population;
vector<double> holdout_fitness_values;
vector<vector< bool > > random_solutions;
vector<double> random_solutions_fitnesses;
vector<vector<int> > all_linkage_sets;
set<vector<int> > set_all_linkage_sets;
int max_all_linkage_pairs_size = 0, max_all_linkage_sets_size = 0;
unordered_set<vector<int>, hash_int_vector> found_linkage_subsets;
double vtr;

unordered_set <vector<int>, hash_int_vector > linkage_pairs;
unordered_map<vector<int>, double, hash_int_vector> wCoef;
Problem *problem_instance;
char instance_path[1000];
int k,s;

void InterpretCommandLine(int argc, char **argv);
void ParseOptions(int argc, char **argv, int *index);
void ParseParameters(int argc, char **argv, int *index);
void InitializeProblem();

void generate_subsets_of_given_size(vector<int>&set, int k, int start, vector<int> &current_combination, vector< vector<int> > &all_combinations);
void generate_all_subsets_of_set(vector<int> &set, vector<vector<int> > &subsets_list);
void calculate_xor(vector<int> &a, vector<int> &b, vector<bool> &solution);
double calculate_probe(vector<vector<int> >&subsets, vector<int> &complement);
void detect_linkage(int N_probes);
void compute_walsh_coefs_after_detect_linkage(vector<vector<int> > &subsets, unordered_map<vector<int>, double, hash_int_vector> &wCoef);
void update_r2();

double approximate(unordered_map<vector<int>, double, hash_int_vector> &coefs);
double get_approximation(int N_probes, bool initialize_linkage_sets, bool calculate_r2);
double get_approximation_heckendorn(int N_probes, bool initialize_linkage_sets, bool calculate_r2);
double perturb_ith_bit(vector<bool> &solution, int i);

void write_statistics(char *fileoption);

double calculate_fitness_change(vector<bool> &prev_solution, vector<bool> &solution_next, vector<int> &bits_changed);
void write_coefficients();
void optimize_surrogate();
bool check_vtr();
void recalc_elitist();
void Run();

void InitializeProblem()
{
  createProblemInstance(problem_index, number_of_variables, &problem_instance, instance_path, k, s);
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

double installedProblemEvaluation(vector <bool> &solution)
{
  double objective_value;

  if (!SKIP_ALREADY_EVALUATED)
  {
    archiveRecord record = checkAlreadyEvaluated(solution, evaluated_solutions);
    bool is_found = record.found;
    double found_value = record.value;    

    if (is_found)
    {
      //printf("FOUND ALREADY EVALUATED! value = %f\n", found_value);
      //number_of_evaluations++;
      return found_value;
    }
  }

  objective_value = problem_instance->calculateFitness(solution);

  if (!SKIP_ALREADY_EVALUATED)
    addToEvaluated(solution, objective_value, evaluated_solutions);
  
  number_of_evaluations++;

  return objective_value;
}


void InterpretCommandLine(int argc, char **argv)
{
  int index = 1;
  ParseOptions(argc, argv, &index);
  ParseParameters(argc, argv, &index);
}

/**
 * Informs the user of an illegal option and exits the program.
 */
void OptionError(char **argv, int index)
{
  printf("Illegal option: %s\n\n", argv[index]);
}

/**
 * Parses only the parameters from the command line.
 */

void ParseOptions(int argc, char **argv, int *index)
{
  double dummy;

  USE_HECKENDORN_ALGORITHM = 0;
  SKIP_ALREADY_EVALUATED = 0;

  for(; (*index) < argc; (*index)++)
  {
    if(argv[*index][0] == '-')
    {
      /* If it is a negative number, the option part is over */
      if(sscanf(argv[*index], "%lf", &dummy) && argv[*index][1] != '\0')
        break;

      if(argv[*index][1] == '\0')
        OptionError(argv, *index);
      else if(argv[*index][2] != '\0')
        OptionError(argv, *index);
      else
      {
        switch(argv[*index][1])
        {
          case 'H': USE_HECKENDORN_ALGORITHM = 1; break;
          case 'N': SKIP_ALREADY_EVALUATED = 1; break;
          default : OptionError(argv, *index);
        }
      }
    }
    else /* Argument is not an option, so option part is over */
     break;
  }

  cout << "USE_HECKENDORN_ALGORITHM:" << USE_HECKENDORN_ALGORITHM << "\n";
  
}

void ParseParameters(int argc, char **argv, int *index)
{
  int noError;
  
  noError = 1;
  int cnt = 0;
  noError = noError && sscanf(argv[*index+cnt], "%d", &problem_index);
  cnt++;
  noError = noError && sscanf(argv[*index+cnt], "%d", &number_of_variables);
  cnt++;
  noError = noError && sscanf(argv[*index+cnt], "%d", &MAX_ORDER);
  cnt++;
  noError = noError && sscanf(argv[*index+cnt], "%lf", &R2_THRESHOLD);
  cnt++;
  noError = noError && sscanf(argv[*index+cnt], "%s", &folder);
  cnt++;
  noError = noError && sscanf(argv[*index+cnt], "%d", &OPTIMIZATION_ALGORITHM);
  cnt++;

  noError = noError && sscanf(argv[*index+cnt], "%d", &timeLimitMilliseconds);
  cnt++;
  noError = noError && sscanf(argv[*index+cnt], "%d", &timeLimitMillisecondsOptimization);
  cnt++;
  noError = noError && sscanf(argv[*index+cnt], "%lld", &RANDOM_SEED);
  cnt++;
  
  if (problem_index == 1 || problem_index == 2)
  {
    noError = noError && sscanf(argv[*index+cnt], "%d", &k);
    cnt++;
    noError = noError && sscanf(argv[*index+cnt], "%d", &s);
    cnt++;
  }

  if (problem_index == 2 || problem_index == 3 || problem_index == 4)
  {
    noError = noError && sscanf(argv[*index+cnt], "%s", &instance_path);
    cnt++;
  }

  //if (MAX_ORDER == -1)
  //  MAX_ORDER = number_of_variables;

  if(!noError)
    cout << "Error parsing parameters." << endl;

  cout << "problem_index:" << problem_index << endl;
  cout << "number_of_variables:" << number_of_variables << endl;
  cout << "MAX_ORDER:" << MAX_ORDER << endl;
  cout << "R2_THRESHOLD:" << R2_THRESHOLD << endl;
  cout << "folder:" << folder << endl;
  cout << "time limit for linkage learning (sec.) :" << timeLimitMilliseconds/1000 << endl;
  cout << "time limit for optimization (sec.) :" << timeLimitMillisecondsOptimization/1000 << endl;
  cout << "OPTIMIZATION_ALGORITHM :" << OPTIMIZATION_ALGORITHM << endl;
  cout << "RANDOM_SEED :" << RANDOM_SEED << endl;

  if (problem_index == 2 || problem_index == 3 || problem_index == 4)
    cout << "instance path:" << instance_path << endl;
  if (problem_index == 1 || problem_index == 2)
    cout << "problem parameters: subfunctions size: " << k << " step: " << s << endl;

}


/**
 * Checks whether the selected options are feasible.
 */
void CheckOptions()
{
  if(number_of_variables < 1)
  {
    printf("\n");
    printf("Error: number of genes < 1 (read: %d). Require number of genes >= 1.", number_of_variables);
    printf("\n\n");

    exit(0);
  }
}

void generate_subsets_of_given_size(vector<int>&set, int k, int start, vector<int> &current_combination, vector< vector<int> > &all_combinations)
{ 
  if (current_combination.size() == k)
  {
    all_combinations.push_back(current_combination);
    return;
  }
  else if (start >= set.size())
    return;

  current_combination.push_back(set[start]);
  generate_subsets_of_given_size(set, k, start+1, current_combination, all_combinations);
  current_combination.pop_back();

  generate_subsets_of_given_size(set, k, start+1, current_combination, all_combinations);
}

void generate_subsets_of_set_of_given_size(vector<int> &set, vector<vector<int> > &subsets_list, int k)
{
  subsets_list.clear();
  vector<int> current_combination;
  generate_subsets_of_given_size(set, k, 0, current_combination, subsets_list);
}

void generate_all_subsets_of_set(vector<int> &set, vector<vector<int> > &subsets_list)
{
    subsets_list.clear();
    for (int cardinality=0; cardinality<=set.size(); ++cardinality)
    {
      vector<int> current_combination;
        generate_subsets_of_given_size(set, cardinality, 0, current_combination, subsets_list);
    }
    //print_vector_of_vectors(subsets_list);
}

void calculate_xor(vector<int> &a, vector<int> &b, vector<bool> &solution)
{
  for(int i = 0; i < a.size(); ++i)
    solution[a[i]] = 1;

  for(int i = 0; i < b.size(); ++i)
    solution[b[i]] = (solution[b[i]] + 1) % 2;
}


double calculate_probe(vector<vector<int> > &mask_subsets, vector<int> &complement)
{
    double probe_value = 0.0;
    vector<bool> tmp_solution(number_of_variables, 0);
        
    for (int i = 0; i < mask_subsets.size(); ++i)
    {
        int mask_order = mask_subsets[i].size();
        int coef;
        if (mask_order % 2 == 0)
          coef = 1;
        else
          coef = -1;

        fill(tmp_solution.begin(), tmp_solution.end(), 0);
        calculate_xor(mask_subsets[i], complement, tmp_solution);
        double f = installedProblemEvaluation(tmp_solution);
        //print_solution(solution);
        //cout << coef << " " << f << endl;
        probe_value += coef * f;
    }

    return probe_value;
}

int find_linkage(vector<bool> &solution1, vector<bool> &solution2, int v)
{
  vector<int> distinct_positions;
  for (int i = 0; i < number_of_variables; ++i)
  {
    if (i == v) continue;
    if (solution1[i] != solution2[i])
      distinct_positions.push_back(i);
  }

  if (distinct_positions.size() == 0)
    return -1;
  if (distinct_positions.size() == 1)
    return distinct_positions[0];
  

  vector<bool> solution3(solution1.begin(), solution1.end());
  for (int i = 0; i < distinct_positions.size() / 2; ++i)
  {
    solution3[distinct_positions[i]] = solution2[distinct_positions[i]];
  }

  double perturb_s1 = perturb_ith_bit(solution1, v);
  double perturb_s3 = perturb_ith_bit(solution3, v);
  //cout << perturb_s1 << " " << perturb_s3 << endl;
  if (fabs(perturb_s1 - perturb_s3) > EPS)
    return find_linkage(solution1, solution3, v);
  else  
    return find_linkage(solution2, solution3, v);
}

double perturb_ith_bit(vector<bool> &solution, int i)
{
  double f1 = installedProblemEvaluation(solution);
  solution[i] = !solution[i];
  double f2 = installedProblemEvaluation(solution);
  solution[i] = !solution[i];
  return f1-f2;  

}

void detect_linkage(int N_probes)
{
  if (GetMilliSecondsRunning(timestamp_start) - gomea_time > timeLimitMilliseconds)
  {
    write_coefficients();
    write_statistics("w");
    optimize_surrogate();
    write_statistics("w");
    recalc_elitist();
    exit(0);
  }

  vector<int> linkage_pair (2);
  vector<bool> s1 (number_of_variables), s2 (number_of_variables);
     
  int valid_tests = 0, tests = 0;
    
  for(int i = 0; i < number_of_variables; ++i)
  {
    int valid_tests = 0;

    //while(valid_tests < N_probes)
    for (int p =0; p < N_probes;++p)
    {      
  
      generate_random_solution(s1);
      generate_random_solution(s2);      

      int deps = 0;
      for (unordered_set <vector<int>, hash_int_vector >::iterator it = linkage_pairs.begin(); it != linkage_pairs.end(); ++it)
      {
        vector<int> pair = *it;
        if (pair[0] == i)
        {
          s2[pair[1]] = s1[pair[1]]; 
          deps++;
        }
        if (pair[1] == i)
        {
          s2[pair[0]] = s1[pair[0]];         
          deps++;
        }
      }
      s2[i] = s1[i];

      bool valid_test = false;
      for (int k = 0; k < number_of_variables; ++k)
      {
        if (s1[k] != s2[k])
          valid_test = true;
      }
      if (!valid_test)
      {
        tests++;
        continue;
      }
      valid_tests++;

      //cout <<i << "perturabation\n";
      double perturb_s1 = perturb_ith_bit(s1, i);
      double perturb_s2 = perturb_ith_bit(s2, i);
      

      if (fabs(perturb_s1 - perturb_s2) > EPS)
      {
        //cout <<"finding linkage started\n";
        int j = find_linkage(s1, s2, i);
        //cout << i << " " << j << endl;
        if (j >= 0)
        {
          linkage_pair[0] = i;
          linkage_pair[1] = j;

          sort(linkage_pair.begin(), linkage_pair.end());
          linkage_pairs.insert(linkage_pair);
        }
        //cout <<"finding linkage finished\n";
      }
    }
  }
  //cout << "valid tests:" << valid_tests << endl;
}


bool check_linkage_set_by_probes(vector<int> &mask, int N_probes)
{
  if (GetMilliSecondsRunning(timestamp_start) - gomea_time > timeLimitMilliseconds)
  {
    write_coefficients();
    optimize_surrogate();
    write_statistics("w");
    recalc_elitist();
    exit(0);
  }

  if (found_linkage_subsets.find(mask) != found_linkage_subsets.end())
    return 1;

  vector<vector<int > > mask_subsets;
  generate_all_subsets_of_set(mask, mask_subsets);

  int min_index = 0, max_index = N_probes;

  vector<bool> boolean_mask(number_of_variables, 0);
  for (int i = 0; i < mask.size(); ++i)
    boolean_mask[mask[i]] = 1;

  for (int i = min_index; i < max_index; ++i)
  {
    //vector<int> probe_to_check(mask.begin(), mask.end());

    vector<bool> solution(number_of_variables);
    
    generate_random_solution(solution);
    
    vector<int> int_solution, mask_complement;
  
    for (int k = 0; k < number_of_variables; ++k)
    {
      if (solution[k]==1 && boolean_mask[k] == 0)
        mask_complement.push_back(k);
    }
    
    double probe_value;
    double multiplier = 1.0 / pow(2.0, mask.size());
     probe_value = multiplier * calculate_probe(mask_subsets, mask_complement);
      
    if (mask_complement.size() == 0)
      all_zeros_probes[mask] = probe_value;
     
    if (fabs(probe_value) > 1e-6)
    { 
     found_linkage_subsets.insert(mask);
     return 1;
   }
  }
  return 0;
}

void compute_walsh_coefs_after_detect_linkage(vector<vector<int> > &all_linkage_sets, unordered_map<vector<int>, double, hash_int_vector> &wCoef)
{
  map<vector<int>, double> hypergraph;

  sort(all_linkage_sets.begin(), all_linkage_sets.end(), vec_comparison_by_size());
  for (int i = 0; i < all_linkage_sets.size(); ++i)
  {
    vector<int> subset = all_linkage_sets[i];
    double multiplier = 1.0 / pow(2, subset.size());
    double all_zeros_probe;

    if (all_zeros_probes.find(subset) != all_zeros_probes.end())
    {
        all_zeros_probe = all_zeros_probes[subset];
    }
    else
    {
      vector<vector<int> > mask_subsets;
      generate_all_subsets_of_set(subset, mask_subsets);      
      vector<int> complement;

      all_zeros_probe = multiplier * calculate_probe(mask_subsets, complement);
      all_zeros_probes.insert(pair<vector<int>, double>(subset, all_zeros_probe));
    }
    hypergraph.insert(pair<vector<int>, double>(subset, all_zeros_probe));
  }

  for (int i = all_linkage_sets.size()-1; i >= 0; --i)
  {
    vector<int> m = all_linkage_sets[i];
    double probeValue = hypergraph[m];

    if (wCoef.find(m) != wCoef.end())
        wCoef[m] = wCoef[m] + probeValue;
    else
        wCoef[m] = probeValue;


    vector<vector<int> > subsets;
    generate_all_subsets_of_set(m, subsets);
      
    for (int j = 0; j < subsets.size(); ++j)
    { 
      vector<int> a = subsets[j];

      if (a.size() == m.size())
          continue;
      
      if (wCoef.find(a) != wCoef.end())
          wCoef[a] = wCoef[a] - wCoef[m];
      else
          wCoef[a] = -wCoef[m];
    }
  }
}


double approximate(unordered_map<vector<int>, double, hash_int_vector> &coefs)
{
    vector<double> all_preds;

    for (int i = 0; i < holdout_population.size(); ++i)
    {
        vector<bool>solution = holdout_population[i];
        double pred = 0.0;
        for (unordered_map<vector<int>, double, hash_int_vector>::iterator it = coefs.begin(); it != coefs.end(); ++it)
        {
            int dot_product = 0;
            for (int j = 0; j < it->first.size(); ++j)
            {
              int ind = it->first[j];
              dot_product += int(solution[ind]);
            }
            if (dot_product % 2 == 0)
              dot_product = 1;
            else
              dot_product = -1;

            pred += it->second * dot_product;
        }
        //cout << holdout_fitness_values[i] << " ~" << pred << endl;
        all_preds.push_back(pred);
    }
    double r2 = r2_score(holdout_fitness_values, all_preds);
    return r2;
}


void update_r2()
{
  //print_vector_of_vectors(all_linkage_sets);
  cout << "NUMBER OF COEFFICIENTS:" << all_linkage_sets.size() << endl;
  cout << "making approximation\n";
  cout << "before compute walsh: " << number_of_evaluations << endl;
  
  wCoef.clear();
  compute_walsh_coefs_after_detect_linkage(all_linkage_sets, wCoef);
  double r2 = approximate(wCoef);    
  
  cout << "updated r2: " << r2 << " "  << endl;
  if (r2 > best_r2)
     best_r2 = r2;
}

double get_approximation(int N_probes, bool initialize_linkage_sets, bool calculate_r2)
{
  //all_linkage_sets.clear();

  if (initialize_linkage_sets)
  {
    vector<int> vec;
    set_all_linkage_sets.insert(vec);
    all_linkage_sets.push_back(vec);
    for (int i = 0; i < number_of_variables; ++i)
    {
      vector<int> vec;
      vec.push_back(i);
      set_all_linkage_sets.insert(vec);
      all_linkage_sets.push_back(vec);
    }
  }

  //print_vector_of_vectors(all_linkage_sets);

  detect_linkage(N_probes);

  vector<vector<vector<int> > >levels;
  levels.resize(number_of_variables+1);
  //cout<<"-----";
  for (unordered_set <vector<int>, hash_int_vector > ::iterator it = linkage_pairs.begin(); it != linkage_pairs.end(); ++it)
  {
    levels[2].push_back(*it);
    if (set_all_linkage_sets.find(*it) == set_all_linkage_sets.end())
    {
      all_linkage_sets.push_back(*it);
      set_all_linkage_sets.insert(*it);
    }
  }

  pair<int, int> pair;

  if (linkage_pairs.size() > max_all_linkage_pairs_size)
  {

    max_all_linkage_pairs_size = linkage_pairs.size();

    int order_limit = MAX_ORDER;
    if (MAX_ORDER == -1)
      order_limit = number_of_variables;
    
    for (int level = 3; level <= order_limit; ++level)
    {
      //cout << level << endl;

      vector<vector<int> > current_level;       
      vector<vector<int> > prev_level = levels[level-1];
      //cout<<"prev level" <<endl;
      //print_vector_of_vectors(prev_level);

      for (int i = 0; i < prev_level.size(); ++i)
      {
        vector<int> current_set = prev_level[i];
        //print_vector(current_set);

        for (int j = current_set[current_set.size()-1] + 1; j < number_of_variables; ++j)
        {
          int flag = 1;
          for (int k = 0; k < current_set.size(); ++k)
          {
              vector<int> pair(2);
              pair[0] = j;
              pair[1] = current_set[k];
              sort(pair.begin(), pair.end());
              //cout<<pair.first << " " << pair.second << endl;
              //cout << int(find(pairs.begin(), pairs.end(), pair) == pairs.end()) << endl;

              if (linkage_pairs.find(pair) == linkage_pairs.end())
              {
                  flag = 0;
                  break;
              }
          }
          if (flag == 0)
              continue;
          //cout << "find " << j << " " <<flag << endl;
          current_set.push_back(j);
          //print_vector(current_set);
          current_level.push_back(current_set);

          if (set_all_linkage_sets.find(current_set) == set_all_linkage_sets.end())
          {
            all_linkage_sets.push_back(current_set);
            set_all_linkage_sets.insert(current_set);
          }

          current_set.pop_back();
        }
      }

      cout << "level:  " << level << " " << current_level.size() << endl;
      //print_vector_of_vectors(current_level);
      levels[level] = current_level;
      if (current_level.size() == 0)
        break;
    }
    
    update_r2();
  }

  if (calculate_r2)
    update_r2();

  return best_r2;
}

double get_approximation_heckendorn(int N_probes, bool initialize_linkage_sets, bool calculate_r2)
{
  all_linkage_sets.clear();

  vector<vector<vector<int> > >levels;
  levels.resize(number_of_variables+1);

  vector<int> vec;
  all_linkage_sets.push_back(vec);
    
  if (initialize_linkage_sets)
    found_linkage_subsets.insert(vec);
  
  for (int i = 0; i < number_of_variables; ++i)
  {
    vector<int> vec;
    vec.push_back(i);
    levels[1].push_back(vec);

    all_linkage_sets.push_back(vec);
      
    if (initialize_linkage_sets)
      found_linkage_subsets.insert(vec);  
  }
  
  int order_limit = MAX_ORDER;
  if (MAX_ORDER == -1)
    order_limit = number_of_variables;
  
  int max_order = 2;
  while (true)
  {      
    vector<vector<int> > current_level;       
    vector<vector<int> > prev_level = levels[max_order-1];
    
    unordered_set<vector<int>, hash_int_vector> check_current_level;
    for (int i = 0; i < prev_level.size(); ++i)
    {
      for (int j = prev_level[i][prev_level[i].size()-1] + 1; j < number_of_variables; ++j)
      {
        vector<vector<int> > subsets;
        prev_level[i].push_back(j);

        //print_vector(prev_level[i]);
        generate_subsets_of_set_of_given_size(prev_level[i], subsets, prev_level[i].size() - 1);
        int flag = 1;
        for (int k = 0; k < subsets.size(); ++k)
        {
          if (subsets[k].size() == prev_level[i].size()) continue;

          if (found_linkage_subsets.find(subsets[k]) == found_linkage_subsets.end())
          {
            //print_vector(subsets[k]);
            flag = 0;
            break;
          }
        }
        //print_vector_of_vectors(subsets);
        if (flag == 1)
          check_current_level.insert(prev_level[i]);
        
        prev_level[i].pop_back();        
      }
    }

    cout << max_order << "to check: " << check_current_level.size() << endl;

    for (unordered_set<vector<int>, hash_int_vector> ::iterator it = check_current_level.begin(); it != check_current_level.end(); ++it)
    {
      vector<int> current_set = *it;
      //print_vector(current_set);        
        
      if (check_linkage_set_by_probes(current_set, N_probes))
      {
        //print_vector(current_set);
        current_level.push_back(current_set);
        all_linkage_sets.push_back(current_set);        
      }
    }
    
    levels[max_order] = current_level;
    if (current_level.size() == 0)
      break;

    max_order++;    
  }

  //print_vector_of_vectors(all_linkage_sets);
  cout << "NUMBER OF COEFFICIENTS:" << all_linkage_sets.size() << endl;

  if (all_linkage_sets.size() > max_all_linkage_sets_size)
  {
    cout << "making approximation\n";

    max_all_linkage_sets_size = all_linkage_sets.size();

    update_r2();
  }
  
  if (calculate_r2)
    update_r2();

  return best_r2;
}



void generate_random_population_with_fitness(int size, vector<vector<bool> > &population, vector<double> &fitness_values)
{
  while (population.size() < size)
  {
    vector<bool> solution(number_of_variables);
    generate_random_solution(solution);
    population.push_back(solution);
    double fitness = installedProblemEvaluation(solution);
    fitness_values.push_back(fitness);
  }
}

void Run(int first_run)
{
  cout << "start run" << endl;
  int N_probes = 1;  int first_probe_to_check = 0;
  best_r2 = -1e+308;
  bool r2_need_to_init = true;
  bool first_approximation = false;
  if (first_run)
    first_approximation = true;
  while(true)
  {
      double r2_holdout; 
      if (!USE_HECKENDORN_ALGORITHM)
        r2_holdout = get_approximation(N_probes, first_approximation, r2_need_to_init);
      else
        r2_holdout = get_approximation_heckendorn(N_probes, first_approximation, r2_need_to_init);

      first_approximation = false;
      r2_need_to_init = false;

      cout << "Population size:" << population.size() << " | N probes:" << N_probes << " | R2 holdout: " << r2_holdout << " " << max_all_linkage_pairs_size << endl;
      total_number_of_probes += N_probes;
      if (r2_holdout >= R2_THRESHOLD)
         break;
  }
}

void write_statistics(char *fileoption)
{
  char       string[10000], full_filename[1000];
  int        i;
  FILE      *file;

  sprintf(full_filename, "%s/linkage_discovery.dat", folder);
  double running_time_milliseconds = GetMilliSecondsRunningSinceTimeStamp(timestamp_start);
  cout << "RUNNING TIME " << running_time_milliseconds / 1000.0 << " SECONDS | EXCLUDING GOMEA:" << (running_time_milliseconds  - gomea_time)/ 1000.0<< endl;
  cout << "EVALUATIONS " << number_of_evaluations  << endl;

  file = fopen(full_filename, fileoption);

  if (fileoption == "w")
  {
    sprintf(string, "Number_of_evaluations full_time,seconds optimization_time,seconds, probes\n", number_of_evaluations, running_time_milliseconds/1000.0, gomea_time/1000.0, total_number_of_probes);
    fputs(string, file);  
  }
  
  sprintf(string, "%ld %.3lf %.3lf %d\n", number_of_evaluations, running_time_milliseconds/1000.0, gomea_time/1000.0, total_number_of_probes);
  fputs(string, file);

  fclose(file);
}

void write_coefficients()
{
  char       string[number_of_variables*10], string1[1000], full_filename[1000];
  int        i;
  FILE      *file;

  sprintf(full_filename, "%s/walsh_coefficients.dat", folder);
  file = fopen(full_filename, "w");
  cout << full_filename << endl;
  for (unordered_map<vector<int>, double, hash_int_vector> :: iterator it = wCoef.begin(); it != wCoef.end(); ++it)
  {
    vector<int> coef = it->first;
    double value = it->second;
    if (fabs(value) < 1e-6)
      continue;

    vector<int> mask(number_of_variables, 0);
    for (int i = 0; i < coef.size(); ++i)
      mask[coef[i]] = 1;
    
    sprintf(string, "");
    int index = 0;
    for (int i = 0; i < number_of_variables; ++i)
    {
      //sprintf(string1, "%d ", mask[i]);
      index += sprintf(&string[index], "%d ", mask[i]);   
      //else
      // strcat(string, "1 ");      
      
    }
    sprintf(&string[index], "%.8lf\n", value);
    fputs(string, file);
  }
    
  fclose(file);
}

void recalc_elitist()
{
  char  c, string_vtr[10000], str_elitist[1000], filename[1000];
  int   k;
  FILE *file;

  vector<bool> solution(number_of_variables, 0);
  sprintf( filename, "%s/elitistsolution.dat",  (char * )folder);
  file = fopen(filename, "r");
  c = fgetc(file);
  k = 0;
  while(c != '\n' && c != EOF)
  {
    if (c == '1')
        solution[k] = 1;
    else
        solution[k] = 0;
    c = fgetc(file);    
    k++;
  }
  for (int i = 0; i < number_of_variables; ++i)
    cout << solution[i];
    
  double new_value = installedProblemEvaluation(solution);
  cout << "\nupdated elitist: " << new_value << endl;
  fclose(file);


  sprintf( filename, "%s/elitist.dat",  (char * )folder);
  file = fopen(filename, "r");
  c = fgetc(file);
  k = 0;
  char evals[1000], milliseconds[1000];
  double value;
  while(c != EOF)
  {
    while(c != ' ')
    {
        evals[k] = (int)c;
        c = fgetc(file);        
        k++;
    }
    evals[k]='\0';
    //printf("EVALS: %s\n", evals);

    c = fgetc(file);
    k = 0;
    while(c != ' ')
    {
        milliseconds[k] = (int)c;
        c = fgetc(file);        
        k++;
    }
    milliseconds[k]='\0';
    //printf("TIME: %s\n", milliseconds);

    while(c != '\n' && c != EOF)
        c = fgetc(file);        
    c = fgetc(file);
    k = 0;
  }
  fclose(file);

  char to_write[10000];
  sprintf(to_write, "%s %s %.6lf", evals, milliseconds, new_value);
  sprintf( filename, "%s/elitist.dat",  (char * )folder);
  file = fopen(filename, "a");
  fputs(to_write, file);
  fclose(file);

  sprintf(to_write, "%.6lf", new_value);
  sprintf( filename, "%s/elitistonly.dat",  (char * )folder);
  file = fopen(filename, "w");
  fputs(to_write, file);
  fclose(file);
  
}

bool check_vtr()
{
  char  c, string_vtr[10000], str_elitist[1000], filename[1000];
  int   i;
  FILE *file;

  /////////////////////////////
  //read elitist
  sprintf( filename, "%s/elitistonly.dat",  (char * )folder);
  file = fopen( filename, "r" );
  if( file == NULL )
    return 1;

  c = fgetc(file);
  k = 0;
  while(c != '\n' && c != EOF)
  {
    str_elitist[k] = (char) c;
    c      = fgetc(file);
    k++;
  }
  str_elitist[k] = '\0';
  double elitist = atof(str_elitist);
  cout << elitist << " " << vtr << " " << fabs(elitist-vtr) << endl;

  if (elitist >= vtr)
    return true;
  return false;
  
}

void optimize_surrogate()
{
  if (timeLimitMillisecondsOptimization >= 3600*24*1000)
    exit(0);

  if (OPTIMIZATION_ALGORITHM == 0)
  {
    char run[1000];
    char walsh_coefficients_filename[1000];
    sprintf(walsh_coefficients_filename, "%s/walsh_coefficients.dat", folder);
    cout << "time limit for GOMEA:" << timeLimitMillisecondsOptimization / 1000.0 << " seconds" << endl;
    sprintf(run, "./BinaryGOMEA -i -r -e -S -v 5 %d 100 %s %d %lld %s", number_of_variables, folder, timeLimitMillisecondsOptimization, RANDOM_SEED, walsh_coefficients_filename);
    printf("%s",run);
    system(run);
  }
  else
  {
    GBOOptimizer *opt = new GBOOptimizer(number_of_variables, OPTIMIZATION_ALGORITHM, folder, wCoef, linkage_pairs);
    try
    {
      opt->optimize(timeLimitMillisecondsOptimization, vtr);
    }
    catch (customException &ex)
    {}

    delete opt;
  }

  timeLimitMillisecondsOptimization *= 2;  
}

int main(int argc, char **argv)
{
  timestamp_start = GetCurrentTimeStampInMilliSeconds();

  long long timestamp_algo;

  InterpretCommandLine(argc, argv);

  srand(RANDOM_SEED);

  InitializeProblem();

  vtr = load_vtr(folder);

  number_of_evaluations = 0;
  
  int holdout_population_size = number_of_variables;
  int k = 2;
  int first_run = 1;
  while (true)
  {    
    cout << holdout_population_size << endl;
  
    generate_random_population_with_fitness(holdout_population_size, holdout_population, holdout_fitness_values);  
    printf("random population for holdout of size %d generated\n", holdout_population_size);

    int skip_run = false;
    if (wCoef.size() > 0)
    {
      double r2 = approximate(wCoef);
      printf("current r2: %.10f\n", r2);
      if (r2 >= R2_THRESHOLD)
        skip_run = true;
    }

    if (skip_run == false)
    {
      Run(first_run);
      first_run = 0;
    }

    holdout_population_size *= 2;
    
    write_statistics("w");
    write_coefficients();

    //optimize_surrogate_with_population_px();
    timestamp_algo = GetCurrentTimeStampInMilliSeconds();
    optimize_surrogate();
    gomea_time += (GetCurrentTimeStampInMilliSeconds() - timestamp_algo);
    cout << "gomea time: " << gomea_time << endl;

    write_statistics("w");
    recalc_elitist();
    if (check_vtr())
      exit(1);
    k++;
  }

  return(0);
}
