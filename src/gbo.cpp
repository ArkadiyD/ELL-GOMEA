#include "gbo.hpp"

GBOOptimizer::GBOOptimizer(int number_of_variables_, int algorithm_, char *folder_, unordered_map<vector<int>, double, hash_int_vector> &wCoef_, unordered_set <vector<int>, hash_int_vector > &linkage_pairs_): 
	number_of_variables(number_of_variables_),
	algorithm(algorithm_),
	folder(folder_),
	wCoef(wCoef_),
	linkage_pairs(linkage_pairs_)
{
  char walsh_coefficients_filename[1000];
  sprintf(walsh_coefficients_filename, "%s/walsh_coefficients.dat", folder.c_str());
  problem_instance = new WalshDecomposition(walsh_coefficients_filename);
  problem_instance->initializeProblem(number_of_variables);

  elitist.genotype.resize(number_of_variables);
  elitist.fitness = -1e+308;


  walsh_coefficients_by_index.resize(number_of_variables);
  for (unordered_map<vector<int>, double, hash_int_vector>::iterator it = wCoef.begin(); it != wCoef.end(); ++it)
  {
    vector<int> subset = it->first;
    double value = it->second;

    walsh_coefficients.push_back(subset);
    walsh_coefficients_values.push_back(value);

    for (int k = 0; k < subset.size(); ++k)
      walsh_coefficients_by_index[subset[k]].push_back(walsh_coefficients.size()-1);
  }
}

void GBOOptimizer::check_elitist_and_vtr(Solution &solution)
{
  /* Update elitist solution */
  if (solution.fitness > elitist.fitness)
  {
  	elitist = solution;
    elitist_number_of_evaluations = number_of_evaluations;
    elitist_time = GetMilliSecondsRunning(timestamp_start);
	  if (solution.fitness >= VTR-1e-6)
    {
        write_elitist_to_file();
        cout << "VTR HIT!" << endl;;
        throw customException("VTR");
    }

    write_elitist_to_file();
  }
}


void GBOOptimizer::write_elitist_to_file()
{
  char       string[10000], filename1[1000], filename2[1000], filename3[1000];
  int        i;
  FILE      *file;

  sprintf(filename1, "%s/elitist.dat", folder.c_str());

  if(!elitist_solution_written_before)
  {
    file = fopen (filename1, "w");
    elitist_solution_written_before = 1;
  }
  else
    file = fopen (filename1, "a");
  /* Write actual values. */
  sprintf( string, "%d %d %.6lf\n", (int)elitist_number_of_evaluations, elitist_time, elitist.fitness);
  fputs( string, file );
  fclose( file );

  sprintf(filename2, "%s/elitistonly.dat", folder.c_str());
  file = fopen( filename2, "w" );
  sprintf( string, "%.6lf", elitist.fitness);
  fputs( string, file );
  fclose( file );

  sprintf(filename3, "%s/elitistsolution.dat", folder.c_str());
  file = fopen( filename3, "w" );
  for (int i = 0; i < number_of_variables; ++i)
  {
  	//cout << elitist[i];

    sprintf( string, "%d", (int)elitist.genotype[i] );
    fputs( string, file );
  }
  //cout << endl;
  fclose( file );

  elitist_solution_written_before = true;
}


double GBOOptimizer::installedProblemEvaluation(Solution &solution)
{
  solution.fitness = problem_instance->calculateFitness(solution.genotype);
  number_of_evaluations += 1;
  check_elitist_and_vtr(solution);

  if (GetMilliSecondsRunning(timestamp_start) > time_limit)
  {
  	cout << "time limit reached!\n";
  	throw customException("time");
  }

  return solution.fitness;
}

double GBOOptimizer::installedProblemEvaluationPartial(Solution &solution, vector<int> &touched_genes_indices, Solution &previous_solution)
{
  double new_evals = 0.0;
  solution.fitness = problem_instance->calculateFitnessPartialEvaluations(solution.genotype, touched_genes_indices, previous_solution.genotype, previous_solution.fitness, new_evals);
  number_of_evaluations += new_evals;

  check_elitist_and_vtr(solution);

  if (GetMilliSecondsRunning(timestamp_start) > time_limit)
  {
  	cout << "time limit reached!\n";
  	throw customException("time");
  }

  return solution.fitness;
}

void GBOOptimizer::dfs(int start_vertex, vector<vector<int> > &graph, vector<int> &labels, int component)
{
  labels[start_vertex] = component;
  for (int i = 0; i < graph[start_vertex].size(); ++i)
  {
    int v = graph[start_vertex][i];
    if (labels[v] == -1)
      dfs(v, graph, labels, component);
  }
}

void GBOOptimizer::find_connected_components(vector<vector<int> > &graph, vector<vector<int> > &components)
{
  vector<int> component_index(number_of_variables, -1);
  int current_component = 0;

  for (int i = 0; i < number_of_variables; ++i)
  {
    if (graph[i].size() == 1 && graph[i][0] == -1) continue;

    if (component_index[i] == -1)
    {
      dfs(i, graph, component_index, current_component);
      //cout << current_component << endl;
      current_component++;
    }
  }

  //print_vector(component_index);

  components.resize(current_component);
  for (int i = 0; i < number_of_variables; ++i)
  {
    if (component_index[i] != -1)
      components[component_index[i]].push_back(i);
  }
  

  //cout << "----------\n";
  //print_vector_of_vectors(components);
  //cout << "----------\n";

}

void GBOOptimizer::partition_crossover(Solution &solution1, Solution &solution2, Solution &child)
{
  vector<vector<int> > graph(number_of_variables);

  for (unordered_set <vector<int>, hash_int_vector > ::iterator it = linkage_pairs.begin(); it != linkage_pairs.end(); ++it)
  {
    vector<int> linkage_pair = *it;
    //print_vector(linkage_pair);
    int v1 = linkage_pair[0];
    int v2 = linkage_pair[1];
    if ((solution1.genotype[v1] != solution2.genotype[v1]) && (solution1.genotype[v2] != solution2.genotype[v2]))
    {
      graph[v1].push_back(v2);
      graph[v2].push_back(v1);
      //cout << "edge: " << v1 << " " << v2 << endl;
    }
  }

  for (int i = 0; i < number_of_variables; ++i)
  {
    if (solution1.genotype[i] == solution2.genotype[i])
    {
      graph[i].resize(1);
      graph[i][0] = -1;     
    }
  }
  // for (int i = 0; i < graph.size(); ++i)
  // {
  //   //cout << i << ":";
  //   //print_vector(graph[i]);
  // }

  vector<vector<int> > components;
  find_connected_components(graph, components);

  child = solution1;
  child.fitness = installedProblemEvaluation(child);

  for (int i = 0; i < components.size(); ++i)
  {
    //double f1 = calculate_fitness_change(solution1, child, components[i]);
    Solution backup = child;

    for (int j = 0; j < components[i].size(); ++j)
    {
      int v = components[i][j];
      if (solution1.genotype[v] == solution2.genotype[v]) continue;
      child.genotype[v] = solution2.genotype[v];
    }
    
    //double f2 = calculate_fitness_change(solution1, child, components[i]);
    child.fitness = installedProblemEvaluationPartial(child, components[i], backup);
    if (backup.fitness > child.fitness || (backup.fitness == child.fitness && random_float_0_1() < 0.5))
      child = backup;

  }

  // /cout << installedProblemEvaluation(solution1) << " " << installedProblemEvaluation(solution2) << " " << installedProblemEvaluation(child) << endl;
}

void GBOOptimizer::hamming_ball_hill_climber(Solution &solution)
{
  //cout << "hmhc starts\n";
  int best_position = -1;
  solution.fitness = installedProblemEvaluation(solution);
  Solution result = solution;
  double best_fitness = result.fitness;

  vector<int> changed(1);

  vector<int> variables_order(number_of_variables);
  iota(variables_order.begin(), variables_order.end(), 0);
  random_shuffle(variables_order.begin(), variables_order.end());

  for (int i = 0; i < number_of_variables; ++i)
  {
    int var = variables_order[i];

    result.genotype[var] = !result.genotype[var];
    changed[0] = var;
    result.fitness = installedProblemEvaluationPartial(result, changed, solution);
    //cout << installedProblemEvaluation(solution) << " " << updated_score << " " << installedProblemEvaluation(result) << endl;

    //print_solution(solution);
    //print_solution(result);

    //assert(installedProblemEvaluation(result)== updated_score);

    if (result.fitness > best_fitness || (result.fitness == best_fitness && random_float_0_1() < 0.5))
    {
      best_fitness = result.fitness;
      best_position = var;
    }

    result.genotype[var] = !result.genotype[var];
  }  

  //cout << start_score << " " << cur_score << " " << best_position << endl;
  if (best_position != -1)
  {
    solution.genotype[best_position] = !solution.genotype[best_position];
    solution.fitness = best_fitness;
  }
}

void GBOOptimizer::perturb(Solution &solution, double alpha)
{
  //print_solution(solution);
  for (int i = 0; i < number_of_variables; ++i)
  {
    if (random_float_0_1() < alpha)
      solution.genotype[i] = !solution.genotype[i];
  }
  //  print_solution(solution);

}

void GBOOptimizer::optimize_surrogate_with_drils(double alpha)
{
  Solution current(number_of_variables);
  generate_random_solution(current.genotype);
  
  hamming_ball_hill_climber(current);

  Solution child(number_of_variables);

  while(true)
  {
    Solution next = current;
    perturb(next, alpha);
    //cout << "next before hc: " << installedProblemEvaluation(next) << endl;
    hamming_ball_hill_climber(next);
    //cout << "next after hc: " << installedProblemEvaluation(next) << endl;

    partition_crossover(current, next, child);
    //cout << "child after PX: " << installedProblemEvaluation(child) << endl;

    if (child.genotype == current.genotype || child.genotype == next.genotype)
    {
      //cout << (int)(child == current) << " " <<  (int)(child == next) << endl;
      //cout << installedProblemEvaluation(next);
      current = next;
    }
    else
    {
      hamming_ball_hill_climber(child);
      current = child;
    }
    //cout << "child after hc: " << installedProblemEvaluation(child) << endl;
  } 
}

void GBOOptimizer::optimize_surrogate_with_hirels()
{
  vector<pair<Solution, int> > stack;
  int evals = 0;

  while(true)
  {
    Solution current(number_of_variables);
    generate_random_solution(current.genotype);
    //cout << "before:" << installedProblemEvaluation(random) << endl;
    hamming_ball_hill_climber(current);
    //cout << "after:" << installedProblemEvaluation(current) << endl;

    int current_level = 0;

    if (stack.size() == 0 or stack[stack.size()-1].second > 0)
      stack.push_back(make_pair(current, current_level));
    else
    {
      //cout << "stack:" << stack.size() << endl;   
      //for (int i = 0; i < stack.size(); ++i)
      //cout << stack[i].second << "|" << stack[i].first.fitness << " ";
      //cout << endl;
      
      //for (int i = 0; i < stack.size(); ++i)
      // cout <<  installedProblemEvaluation(stack[i].first) << " ";     
      //cout << endl;

      bool pxSucess = true;
      while (stack.size() != 0 && pxSucess && stack[stack.size()-1].second == current_level)
      {
        Solution top = stack[stack.size()-1].first, child(number_of_variables);
        stack.pop_back();

        //cout << current_level << " before px:" << installedProblemEvaluation(current) << " " << installedProblemEvaluation(top) << endl;   
        partition_crossover(top, current, child);
        //cout << "after px:" << installedProblemEvaluation(child) << endl;

        pxSucess = (child.genotype != top.genotype) && (child.genotype != current.genotype);
        if (pxSucess)
        {
          hamming_ball_hill_climber(child);
          current = child;
          //cout << "after hbhc:" << installedProblemEvaluation(current) << endl;
          current_level++;
        }
      }

      if (pxSucess)
        stack.push_back(make_pair(current, current_level));

      //cout << " loop ended\n";
    }


    //cout << stack.size() << endl;
    //double d = installedProblemEvaluation(current);
  } 
}

// void GBOOptimizer::optimize_surrogate_with_population_px()
// {
//   int evals = 0;
//   vector<vector< bool> > population(1000);
//   for (int i = 0; i < population.size(); ++i)
//   {
//     population[i].resize(number_of_variables);
//     generate_random_solution(population[i]);
//   }
    
//   while(true)
//   {
//     vector<vector< bool> > population_offsprings(population.begin(), population.end());
//     for (int i = 0; i < population.size(); ++i)
//     {
//       vector<bool> parent1 = population[rand() % population.size()];
//       vector<bool> parent2 = population[rand() % population.size()];
//       //cout << installedProblemEvaluation(parent1) << " " << installedProblemEvaluation(parent2) << endl;

//       partition_crossover(parent1, parent2, population_offsprings[i]);
//       hamming_ball_hill_climber(population_offsprings[i], population_offsprings[i]);
//       double d = installedProblemEvaluation(population_offsprings[i]);
//       //cout << installedProblemEvaluation(population_offsprings[i]) << endl;
//     }

//     population = population_offsprings;
//   }
// }

void GBOOptimizer::optimize(long long time_limit_, double VTR_)
{
	time_limit = time_limit_;
	VTR = VTR_;
	timestamp_start = GetCurrentTimeStampInMilliSeconds();
  cout << " OPTIMIZATION ALGORTIHM " << algorithm << endl;

	if (algorithm == 1)
		optimize_surrogate_with_hirels();
	 else if (algorithm == 2)
	 	optimize_surrogate_with_drils(0.01);
   else if (algorithm == 3)
    optimize_surrogate_with_drils(0.1);
   
	// else if (algorithm == 3)
	// 	optimize_surrogate_with_population_px();

}