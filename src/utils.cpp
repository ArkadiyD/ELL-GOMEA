#include "utils.hpp"

//////////////////////////////////////////////////////////////////////////////

long GetMilliSecondsRunning(long timestamp_start)
{
  return(GetMilliSecondsRunningSinceTimeStamp(timestamp_start));
}

long GetMilliSecondsRunningSinceTimeStamp(long timestamp)
{
  long timestamp_now, difference;

  timestamp_now = GetCurrentTimeStampInMilliSeconds();

  difference = timestamp_now-timestamp;

  return(difference);
}

long GetCurrentTimeStampInMilliSeconds()
{
  struct timeval tv;
  long   result;

  gettimeofday(&tv, NULL);
  result = (tv.tv_sec * 1000) + (tv.tv_usec / 1000);

  return(result);
}

//////////////////////////////////////////////////////////////////////////////

void print_solution(vector<bool> &vec)
{
  for (int i = 0; i < vec.size(); ++i)
    cout << vec[i] << " ";
  cout << endl;
}

void print_vector(vector<int> &vec)
{
  for (int i = 0; i < vec.size(); ++i)
    cout << vec[i] << " ";
  cout << endl;
}

void print_vector_of_vectors(vector<vector<int> > &subsets_list)
{
  for (int i = 0; i < subsets_list.size(); ++i)
    print_vector(subsets_list[i]);
}

void print_population(vector<vector<bool> > &subsets_list)
{
  for (int i = 0; i < subsets_list.size(); ++i)
  {
    for (int j =0; j < subsets_list[i].size(); ++j)
      cout << subsets_list[i][j] << " ";
    cout<< "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

double mean(vector<double> &vec)
{
  double sum = 0.0;
  for (int i = 0; i < vec.size(); ++i)
    sum += vec[i];
  return sum / vec.size();
}

double var(vector<double> &vec)
{
  double mean_value = mean(vec);
  double sum = 0.0;
  for (int i = 0; i < vec.size(); ++i)
    sum += (vec[i] - mean_value) * (vec[i] - mean_value);
  return sum / vec.size();
}

double r2_score(vector<double> &target, vector<double> &pred)
{
  double sse = 0;
  for (int i = 0; i < target.size(); i++)
    sse += (target[i] - pred[i]) * (target[i] - pred[i]);
  sse /= target.size();
  double var_value = var(target);
  return 1.0 - sse / var_value;
}

//////////////////////////////////////////////////////////////////////////////

double random_float_0_1()
{
  return rand() / (double)RAND_MAX;
}

bool randomBool()
{
   return rand() > (RAND_MAX / 2);
}

void generate_random_solution(vector<bool> &solution)
{
  for (int j = 0; j < solution.size(); ++j)
    solution[j] = randomBool();
}


//////////////////////////////////////////////////////////////////////////////

void addToEvaluated(vector <bool> &solution, double objective_value, unordered_map< vector<bool>, double, hash_vector > &evaluated_solutions)
{
  //printf("ADDING SOLUTION with fitness %f\n", objective_value);
  evaluated_solutions.insert(pair<vector<bool>, double> (solution, objective_value));
}

//check whether the solution is already evalauted, if yes return corresponding fitness values
archiveRecord checkAlreadyEvaluated(vector<bool> &solution, unordered_map< vector<bool>, double, hash_vector > &evaluated_solutions)
{
  archiveRecord res;
  res.found = false;

  unordered_map<vector<bool>, double, hash_vector >::iterator it = evaluated_solutions.find(solution);
  if (it != evaluated_solutions.end())
  {
    res.found = true;
    res.value = it->second;
  }
  
  return res;
}

////////////////////////////////////////////////////////////
double load_vtr(char *folder)
{
  char  c, string_vtr[10000], str_elitist[1000], filename[1000];
  int   i;
  FILE *file;

  /////////////////////////////
  //read vtr
  sprintf( filename, "%s/vtr.txt",  (char *)folder);
  file = fopen( filename, "r" );
  if( file == NULL )
  {
    cout << "VTR file not found!\n";
    exit(1);
  }

  c = fgetc(file);
  int k = 0;
  while(c != '\n' && c != EOF)
  {
    string_vtr[k] = (char) c;
    c      = fgetc(file);
    k++;
  }
  string_vtr[k] = '\0';
  double vtr = atof(string_vtr);
  cout << vtr << endl;
  return vtr;
}