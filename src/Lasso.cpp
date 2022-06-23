#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include <functional>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <string>
#include <unistd.h>
#include <signal.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
using namespace std;


#include <mlpack/core.hpp>
#include <mlpack/methods/lars/lars.hpp>
#include <mlpack/core/cv/k_fold_cv.hpp>
#include <mlpack/core/cv/metrics/mse.hpp>

using namespace mlpack;
using namespace mlpack::regression;
using namespace mlpack::cv;

#include "utils.hpp"
#include "problems.hpp"

double EPS = 1e-6;
double R2_THRESHOLD;
long long RANDOM_SEED;
double vtr;

int number_of_evaluations = 0,
  number_of_variables,
  problem_index,
  MAX_ORDER,
  max_evaluated_solutions = 1000000;

long        timestamp_start,                                                    /* The time stamp in milliseconds for when the program was started. */
            timestamp_finish;                                         /* The time stamp in milliseconds for when the algorithm was started (after problem initialization). */
int64_t random_seed,
    random_seed_changing;
char folder[1000];
long long gomea_time = 0;
long int timeLimitMilliseconds, timeLimitMillisecondsOptimization;

unordered_map< vector<bool>, double, hash_vector > evaluated_solutions;

vector<vector<bool> > population;
vector<double> population_fitness_values;
vector<vector<bool> > holdout_population;
vector<double> holdout_fitness_values;
map<vector<int>, double> wCoef;
Problem *problem_instance;
char instance_path[1000];
int k,s;

arma::mat arma_population;
arma::rowvec arma_population_fitness, arma_holdout_population_fitness;

void InterpretCommandLine(int argc, char **argv);
void ParseOptions(int argc, char **argv, int *index);
void ParseParameters(int argc, char **argv, int *index);

void InitializeRandomNumberGenerator();
void InitializeProblem();

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

void write_coefficients();
void optimize_surrogate();
void write_statistics(char *fileoption);
void finish();

void InitializeProblem()
{
  createProblemInstance(problem_index, number_of_variables, &problem_instance, instance_path, k, s);
}


void catchAlarm(int sig)
{
    cerr << "Time limit exceeded!\n";
    finish();
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


double installedProblemEvaluation(vector <bool> &solution)
{
  double objective_value;

  archiveRecord record = checkAlreadyEvaluated(solution, evaluated_solutions);
  bool is_found = record.found;
  double found_value = record.value;    

  if (is_found)
  {
    //printf("FOUND ALREADY EVALUATED! value = %f\n", found_value);
    //number_of_evaluations++;
    return found_value;
  }

  objective_value = problem_instance->calculateFitness(solution);

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

void PrintUsage()
{
  printf("Usage: BinaryGOMEA [-?] [-P] [-e] [-f] [-i] [-p] [-r] [-v] pro dim vtr mng\n");
  printf("   -?: Prints out this usage information.\n");
  printf("   -P: Prints out a list of all installed optimization problems.\n");
  printf("   -p: Use population instead of random solutions\").\n");
  printf("   -l: Use lasso to estimate walsh coefficients (exact calculation by default)\n");
  
  printf("\n");
  printf("  pro: Index of optimization problem to be solved (minimization).\n");
  printf("  dim: Number of variables.\n");
  printf("  max order: Maximum order of walsh coefficients, -1 means try all possible orders until a perfect approximation is obtained\n");
  printf("  r2 threshold to reach\n");
  printf("  folder: folder of filename to save results.\n");
  
  
  exit(0);
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
          case '?': PrintUsage(); break;
          default : OptionError(argv, *index);
        }
      }
    }
    else /* Argument is not an option, so option part is over */
     break;
  } 
}

void ParseParameters(int argc, char **argv, int *index)
{
  int noError;
  int cnt = 0;

  noError = 1;
  
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


  if (MAX_ORDER == -1)
    MAX_ORDER = number_of_variables;


  if(!noError)
    cout << "Error parsing parameters." << endl;

  cout << "problem_index:" << problem_index << endl;
  cout << "number_of_variables:" << number_of_variables << endl;
  cout << "MAX_ORDER:" << MAX_ORDER << endl;
  cout << "R2_THRESHOLD:" << R2_THRESHOLD << endl;
  cout << "folder:" << folder << endl;
  cout << "time limit for linkage learning (sec.) :" << timeLimitMilliseconds/1000 << endl;
  cout << "time limit for optimization (sec.) :" << timeLimitMillisecondsOptimization/1000 << endl;
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
  //for (int i = 0; i < current_combination.size(); ++i)
  //  cout << current_combination[i] << " ";
  //cout << endl;
  
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

double approximate(map<vector<int>, double> &coefs)
{
    vector<double> all_preds;

    for (int i = 0; i < holdout_population.size(); ++i)
    {
        vector<bool>solution = holdout_population[i];
        double pred = 0.0;
        for (map<vector<int>, double>::iterator it = coefs.begin(); it != coefs.end(); ++it)
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

void optimize_surrogate()
{
  char run[1000];
  char walsh_coefficients_filename[1000];
  sprintf(walsh_coefficients_filename, "%s/walsh_coefficients.dat", folder);
  cout << "time limit for GOMEA:" << timeLimitMillisecondsOptimization / 1000.0 << " seconds" << endl;
  sprintf(run, "./BinaryGOMEA -i -r -e -S -v 5 %d 100 %s %d %lld %s", number_of_variables, folder, timeLimitMillisecondsOptimization, RANDOM_SEED, walsh_coefficients_filename);
  printf("%s",run);
  system(run);

  timeLimitMillisecondsOptimization *= 2;
}


void calculate_walsh_features(vector<vector<bool> > &population, vector<vector<int> > &subsets_list, arma::Mat<double> &matrix)
{
    for (int i = 0; i < subsets_list.size(); ++i)
    {
        for (int k = 0; k < population.size(); ++k)
        {
            int value = 0;

            for (int j = 0; j < subsets_list[i].size(); ++j)
            {
                int ind = subsets_list[i][j];
                if (population[k][ind] == 1)
                    value++;
            }

            if (value % 2 == 0)
                value = 1;
            else
                value = -1;

            matrix(k, i) = value;
        }
    } 
    matrix = matrix.t();  
}

double fit_model(arma::mat &features, arma::rowvec &targets, arma::mat &holdout_features, arma::rowvec &holdout_targets, vector<vector<int> > &subsets_list)
{
    if (GetMilliSecondsRunning(timestamp_start) - gomea_time > timeLimitMilliseconds)
        finish();

    arma::rowvec normed_targets(targets.size());
    //cout << targets << endl;
    double max = -1e+308;
    for (int i = 0; i < targets.size(); ++i)
    {
        if (fabs(targets(i)) > max)
            max = fabs(targets(i));
    }
    for (int i = 0; i < targets.size(); ++i)
        normed_targets(i) = targets(i) / max;
    //cout << "MAX " << max << endl;
    double lambda = 1e-6, best_lambda, min_error=1e+308;

    cout << "start fitting" << endl;
    
    signal(SIGALRM, catchAlarm);
    alarm(timeLimitMilliseconds / 1000); 

    while (lambda <= 1)
    {   
        if (GetMilliSecondsRunning(timestamp_start) - gomea_time > timeLimitMilliseconds)
          finish();

        mlpack::cv::KFoldCV<mlpack::regression::LARS, mlpack::cv::MSE> cv(3, features, normed_targets, true);
        double mse = cv.Evaluate(true, lambda, 1e-6, 1e-8);
        if (mse < min_error)
        {
            min_error = mse;
            best_lambda = lambda;
        }
        //cout << mse << endl;
        lambda *= 10;        
    }

    cout << "cv done" << endl;
    LARS lars(true, best_lambda, 1e-6, 1e-8);
    lars.Train(features, normed_targets);
    arma::rowvec predictions(holdout_targets.size());
    lars.Predict(holdout_features, predictions);
    cout << best_lambda << " " << min_error << endl;
    //cout << predictions;
    //cout << holdout_targets;

    wCoef.clear();
    //cout << "n coefficients" << lars.Beta().size() << endl;
    for (int i = 0; i < subsets_list.size(); ++i)
    {
        //for (int j = 0; j < subsets_list[i].size();++j)
        //    cout << subsets_list[i][j] << " ";
        //cout << lars.Beta()(i)*max << endl;
        wCoef[subsets_list[i]] = lars.Beta()(i) * max;
    }
    
    double r2 = approximate(wCoef);
    cout << "R2:" << r2 << endl;    
    
    return r2;
}

void finish()
{
    write_coefficients();
    write_statistics("w");
    optimize_surrogate();
    write_statistics("w");
    recalc_elitist();
    exit(0);
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
    sprintf(string, "Number_of_evaluations full_time,seconds optimization_time,seconds\n", number_of_evaluations, running_time_milliseconds/1000.0, gomea_time/1000.0);
    fputs(string, file);  
  }
  
  sprintf(string, "%ld %.3lf %.3lf\n", number_of_evaluations, running_time_milliseconds/1000.0, gomea_time/1000.0);
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
  for (map<vector<int>, double> :: iterator it = wCoef.begin(); it != wCoef.end(); ++it)
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


void Run(int first_run)
{
  cout << "start run" << endl;

  if (first_run)
    generate_random_population_with_fitness(number_of_variables, population, population_fitness_values);
  cout << "random population generated\n";

  arma_population_fitness = arma::conv_to<arma::rowvec>::from(population_fitness_values);
  arma_holdout_population_fitness = arma::conv_to<arma::rowvec>::from(holdout_fitness_values);

  vector<int> variables;
  for (int i = 0; i < number_of_variables; ++i)
    variables.push_back(i);

  while(true)
  {
    vector<vector<int> > subsets_list;
    generate_subsets_of_set_of_given_size(variables, subsets_list, 0);
            
    for (int max_order = 1; max_order <= MAX_ORDER; ++max_order)
    {
        generate_subsets_of_set_of_given_size(variables, subsets_list, max_order);
        cout << max_order << " " << subsets_list.size() << " " << population.size() << endl;
        if (subsets_list.size() >= population.size())
            break;

        cout << population.size() << " " << subsets_list.size() << endl;
        //exit(1);
        arma::mat features((int)population.size(), (int)subsets_list.size(), arma::fill::zeros);
        arma::mat holdout_features((int)holdout_population.size(), (int)subsets_list.size(), arma::fill::zeros);
        
        calculate_walsh_features(population, subsets_list, features);
        calculate_walsh_features(holdout_population, subsets_list, holdout_features);
        
        double r2 = fit_model(features, arma_population_fitness, holdout_features, arma_holdout_population_fitness, subsets_list);
        if (r2 >= R2_THRESHOLD)
            return;
    }
    
    generate_random_population_with_fitness(population.size()*2, population, population_fitness_values);
    arma_population_fitness = arma::conv_to<arma::rowvec>::from(population_fitness_values);
  }
}

int calculate_population_size(int L, int k, int M)
{
  double eps = 1e-3;
  double two_power_k = pow(2, k);
  double numerator = log2(1-pow(1-eps, 1.0/(M*two_power_k)));
  double denominator = log2(1-1.0/two_power_k);
  return ceil(numerator / denominator); 
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

    timestamp_algo = GetCurrentTimeStampInMilliSeconds();
    optimize_surrogate();
    gomea_time += (GetCurrentTimeStampInMilliSeconds() - timestamp_algo);
    cout << "gomea time: " << gomea_time << endl;

    write_statistics("w");
    recalc_elitist();
    if (check_vtr())
      exit(0);

  }
 
  return(0);
}
