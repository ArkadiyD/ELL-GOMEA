#pragma once

#include <vector>
#include <sys/time.h>
#include <cstddef>
#include <iostream>
#include <unordered_map>
using namespace std;

struct archiveRecord
{
  bool found;
  double value;
};

struct hash_vector { 
  size_t operator()(const vector<bool> &vec) const
    { 
        hash<vector<bool> > hash_vector_bool; 
      size_t hash = hash_vector_bool(vec); 
        return hash; 
    } 
}; 

struct hash_int_vector { 
  size_t operator()(const vector<int> &vec) const
    { 
      size_t hash_value = 0;
        hash<int> hash_int;
        for(int i = 0; i < vec.size(); ++i)
          hash_value ^= hash_int(vec[i]);

        return hash_value; 
    } 
}; 

struct hash_pair { 
  size_t operator()(const pair<int, int> &pair_) const
    { 
        hash<int> hash_int; 
      size_t hash1 = hash_int(pair_.first);
      size_t hash2 = hash_int(pair_.second);
        return hash1 + hash2; 
    } 
}; 

struct vec_comparison_by_size
{
    inline bool operator() (const vector<int>& vec1, const vector<int>& vec2)
    {
        return (vec1.size() < vec2.size());
    }
};

class customException: public exception
{
private:
    string message;

public:
    customException(string message_) : message(message_) { }
    const char * what () const throw ()
    {
        return message.c_str();
    }
};

long GetMilliSecondsRunning(long timestamp_start);
long GetMilliSecondsRunningAfterInit();
long GetMilliSecondsRunningSinceTimeStamp(long timestamp);
long GetCurrentTimeStampInMilliSeconds();

void print_solution(vector<bool> &vec);
void print_vector(vector<int> &vec);
void print_vector_of_vectors(vector<vector<int> > &subsets_list);
void print_population(vector<vector<bool> > &subsets_list);

double mean(vector<double> &vec);
double var(vector<double> &vec);
double r2_score(vector<double> &target, vector<double> &pred);

bool randomBool();
double random_float_0_1();
void generate_random_solution(vector<bool> &solution);

void addToEvaluated(vector <bool> &solution, double objective_value, unordered_map< vector<bool>, double, hash_vector > &evaluated_solutions);
archiveRecord checkAlreadyEvaluated(vector<bool> &solution, unordered_map< vector<bool>, double, hash_vector > &evaluated_solutions);

double load_vtr(char *folder);