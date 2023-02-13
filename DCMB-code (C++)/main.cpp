//Xianjie Guo, Kui Yu
//10/10/2019
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////includes///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <vector>
#include <list>
#include <queue>
#include <map>
#include <algorithm>

#include "./H/chisq.h"

using std::list;
using std::map;
using std::queue;
using std::vector;

///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////defines//////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#define MAX_LINE_LENGTH 1000000
#define MAX_N_STATES 105
#define MAX_STATE_LENGTH 100

#define N_CASES_PER_DF 5

typedef map<std::string, double> MAP_STRING_DOUBLE;
typedef map<double, int> MAP_DOUBLE_INT;

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////global definitions/////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int n_cases, n_vars, Search_Count, Search_Count2, pc_superset_cache_count; //number of network nodes and sample cases

int max_max_cond_size, k_conditon, label_index; //K_conditon is the maximum test set for conditional independence tests

double K_or, K_and;

double k_tradeoff, alpha, average_perform_time, average_n_G2; //average_n_G2 is the average number of conditional independence tests each node takes

clock_t clock_init;

int *n_states;

const char *data_file_path, *mb_write_in_file_path;

int **data_processed;

int ***sep_cache;

MAP_STRING_DOUBLE *Cond_Size;

vector<vector<int>> pc_superset_cache, or_pc_cache, and_pc_cache;

bool *pc_superset_computed;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////functions' prototypes///////////////////////////
///////////////////////////////////////////////////////////////////////////////

void DCMB(int target);

///////////////////////////////////////////////////////////////////////////////

void parse_data(void);
int state_index(int var, char *state, char ***states);
void pc_superset(int target, int *pc, int **sep, vector<int> &or_pc, vector<int> &and_pc);
double min_dep(int var, int target, int compulsory, int *pc, int *sep);
double min_dep_fast(int var, int target, int compulsory, int *pc, int *sep);
int next_cond_index(int n_pc, int cond_size, int *cond_index);
double compute_dep(int var, int target, int *cond);
double Search_dep(int var, int target, int *cond);
int compute_rows(const char *path);
int compute_columns(const char *path);
void report_mb(int target, int *mb);

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////functions' body////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////

int main(int argc, char *argv[])
{
    int i, j;
    const char *algorithm;
    vector<int> targets;
    std::fstream _file;

    if (argc!=3) {
        std::cout<<"Please enter three parameters！"<<std::endl;
        exit(0);
    }

    /* Parameter input */
    data_file_path =argv[1]; // data_file_path = "data/Alarm_s500.txt";
    label_index =argv[2]; // label column, label_index in [1, #variables];

    alpha = 0.01;
    algorithm = "DCMB";
    K_or=0.3;
    K_and=0.3;
    k_conditon = 3;

    // alpha=atof(argv[3]); // alpha = 0.01;
    // algorithm=argv[4]; // algorithm = "DCMB";
    // K_or=atof(argv[5]); // K_or=0.5;
    // K_and=atof(argv[6]); // K_and=0.5;
    // k_conditon=atoi(argv[7]); // k_conditon = 3;

    n_vars = compute_columns(data_file_path);
    n_cases = compute_rows(data_file_path);

    std::cout << "algorithm:" << algorithm << "  "
              << "alpha:" << alpha << "  "
              << "n_cases:" << n_cases << "  "
              << "n_vars:" << n_vars << " "
              << "K_or:" << K_or << " "
              << "K_and:" << K_and << " " << std::endl;

    // if (!strcmp(argv[5],"all")||!strcmp(argv[5],"")) {
    //     for (int var = 0; var < n_vars; ++var) {
    //          targets.push_back(var);
    //     }
    // }
    // else
    // {
    //     const char *sep = ",";
    //     char *p;
    //     p = strtok(argv[5], sep);
    //        while(p){
    //            targets.push_back(atoi(p));
    //            p = strtok(NULL, sep);
    //        }
    // }

    // for (i = 0; i < n_vars; ++i)
    // {
    //     targets.push_back(i);
    // }

    targets.push_back(label_index-1);
    

    Search_Count = 0;
    Search_Count2 = 0;
    pc_superset_cache_count = 0;

    pc_superset_cache.resize(n_vars);
    or_pc_cache.resize(n_vars);
    and_pc_cache.resize(n_vars);
    pc_superset_computed = new bool[n_vars];

    for (i = 0; i < n_vars; ++i)
    {
        pc_superset_computed[i] = false;
    }

    srand((unsigned)time(NULL));

    mb_write_in_file_path = "mb.out";

    average_perform_time = (double)0.0;
    average_n_G2 = (double)0.0;

    k_tradeoff = 1.0;

    n_states = new int[n_vars];

    data_processed = new int *[n_cases];
    for (i = 0; i < n_cases; i++)
        data_processed[i] = new int[n_vars];

    parse_data();

    sep_cache = new int **[n_vars];
    for (i = 0; i < n_vars; ++i)
    {
        sep_cache[i] = new int *[n_vars];
        for (j = 0; j < n_vars; ++j)
        {
            sep_cache[i][j] = new int[max_max_cond_size];
        }
    }

    Cond_Size = new MAP_STRING_DOUBLE[max_max_cond_size + 1];
    for (i = 0; i < max_max_cond_size + 1; i++)
    {
        Cond_Size[i].clear();
    }

    //select algorithm
    if (!strcmp(algorithm, "DCMB"))
    {
        for (unsigned long j = 0; j < targets.size(); ++j)
        {
            clock_init = clock();
            DCMB(targets[j]);
        }
    }
    else
    {
        std::cout << "Error: invalid algorithm name" << std::endl;
        fflush(stdout);
        exit(0);
    }

    average_perform_time = average_perform_time / (double)targets.size();
    average_n_G2 = average_n_G2 / (double)targets.size();


    _file.open("time.out", std::ios::out | std::ios::app); //read from memory into disk and open file in append mode

    if (!_file)
    {
        std::cout << "error" << std::endl;
        return 0;
    }

    _file << average_perform_time << std::endl;
    _file << std::endl;
    _file.close();

    std::cout << "time:" << average_perform_time << std::endl;

    targets.clear();
    return 0;
}

/////////////////////////////////////

void parse_data(void)
{
    int i, j, min_n_states;
    char buffer[(int)MAX_LINE_LENGTH];
    char *token;
    const char *separators = " ,\t\n\r";
    FILE *f_in;
    char ***states; //Two-dimensional string type

    states = new char **[n_vars];
    for (i = 0; i < n_vars; i++)
    {
        states[i] = new char *[(int)MAX_N_STATES];
        for (j = 0; j < (int)MAX_N_STATES; j++)
            states[i][j] = new char[(int)MAX_STATE_LENGTH];
    }

    for (i = 0; i < n_vars; i++)
        n_states[i] = 0;

    f_in = fopen(data_file_path, "r");

    if (f_in == NULL)
    {
        std::cout << "data file open error!" << std::endl;
        exit(0);
    }

    for (i = 0; i < n_cases; i++)
    {
        fgets(buffer, (int)MAX_LINE_LENGTH, f_in);
        if ((int)strlen(buffer) == (int)MAX_LINE_LENGTH - 1)
        {
            printf("\n\n Exceeding MAX_LINE_LENGTH. \n\n");
            exit(0);
        }

        token = strtok(buffer, separators);
        for (j = 0; j < n_vars; j++)
        {
            if ((int)strlen(token) > (int)MAX_STATE_LENGTH - 1)
            {
                printf("\n\n Exceeding MAX_STATE_LENGTH. \n\n");
                exit(0);
            }

            data_processed[i][j] = state_index(j, token, states);
            //**data is not the data in the text, but the data from 0 after processing

            token = strtok(0, separators); //0 can also be changed to NULL to split the remaining strings.
        }
    }

    fclose(f_in);

    min_n_states = 0;
    for (i = 0; i < n_vars; i++)
        if (n_states[i] > 1 && (min_n_states == 0 || n_states[i] < min_n_states))
            min_n_states = n_states[i];

    if (min_n_states == 0)
    {
        printf("\n\n Only one state for all the variables. \n\n");
        exit(0);
    }

    max_max_cond_size = (int)(log((double)n_cases / (double)((int)N_CASES_PER_DF * (min_n_states - 1) * (min_n_states - 1))) / log((double)min_n_states)); //min_n_states may not be representative.
    //If the condition set of the independence test exceeds max_max_cond_size, the test result is no longer reliable.
    for (i = 0; i < n_vars; i++)
    {
        for (j = 0; j < (int)MAX_N_STATES; j++)
            delete[] states[i][j];
        delete[] states[i];
    }
    delete[] states;
}

/////////////////////////////////////

int state_index(int var, char *state, char ***states)
{
    int i;

    for (i = 0; i < n_states[var]; i++)
        if (strcmp(states[var][i], state) == 0)
            return i;

    if (i > (int)MAX_N_STATES - 1)
    {
        printf("\n\n Exceeding MAX_N_STATES. \n\n");
        exit(0);
    }

    strcpy(states[var][i], state);
    n_states[var]++;

    return i;
}

/////////////////////////////////////


/////////////////////////////////////

void DCMB(int target)
{
    int i, j, k, in_cond, n_conds;
    int *pc, *pc2, *pc3, *mb, *cond;
    int **sep, **sep2;
    vector<int> or_pc, and_pc, or_pc2, and_pc2;
    double aux;

    pc = new int[n_vars];
    pc2 = new int[n_vars];
    pc3 = new int[n_vars];
    mb = new int[n_vars];

    cond = new int[max_max_cond_size];

    sep = new int *[n_vars];
    sep2 = new int *[n_vars];
    for (i = 0; i < n_vars; i++)
    {
        sep[i] = new int[max_max_cond_size];
        sep2[i] = new int[max_max_cond_size];
    }

    pc_superset(target, pc, sep, or_pc, and_pc);

    for (vector<int>::iterator iter = or_pc.begin(); iter < or_pc.end(); iter++)
    {
        pc_superset(*iter, pc2, sep2, or_pc2, and_pc2);
        if (pc2[target] == 1)
        {
            pc[*iter] = 1;
        }
    }
    for (vector<int>::iterator iter = and_pc.begin(); iter < and_pc.end(); iter++)
    {
        pc_superset(*iter, pc2, sep2, or_pc2, and_pc2);
        if (pc2[target] == 0)
        {
            pc[*iter] = 0;
        }
        for (i = 0; i < max_max_cond_size; i++)
        {
            sep[*iter][i] = sep2[target][i];
        }
    }

    for (i = 0; i < n_vars; i++)
        mb[i] = pc[i];

    for (i = 0; i < n_vars; i++)
        if (pc[i] == 1)
        {
            pc_superset(i, pc2, sep2, or_pc, and_pc);

            for (vector<int>::iterator iter = or_pc.begin(); iter < or_pc.end(); iter++)
            {
                pc_superset(*iter, pc3, sep2, or_pc2, and_pc2);
                if (pc3[i] == 1)
                {
                    pc2[*iter] = 1;
                }
            }
            for (vector<int>::iterator iter = and_pc.begin(); iter < and_pc.end(); iter++)
            {
                pc_superset(*iter, pc3, sep2, or_pc2, and_pc2);
                if (pc3[i] == 0)
                {
                    pc2[*iter] = 0;
                }
            }

            for (j = 0; j < n_vars; j++)
                if (pc2[j] == 1 && pc[j] == 0 && j != target && mb[j] == 0)
                {
                    in_cond = 0;
                    n_conds = 0;
                    for (k = 0; k < max_max_cond_size; k++)
                    {
                        cond[k] = sep[j][k];

                        if (cond[k] == i)
                            in_cond = 1;

                        if (cond[k] != -1)
                            n_conds++;
                    }

                    if (in_cond == 0 && n_conds < max_max_cond_size)
                    {
                        cond[n_conds] = i;

                        aux = Search_dep(j, target, cond);
                        if (aux >= (double)1.0)
                        {
                            mb[j] = 1;
                        }
                    }
                }
        }

    report_mb(target, mb);

    delete[] mb;
    delete[] pc;
    delete[] pc2;
    delete[] pc3;
    delete[] cond;
    for (i = 0; i < n_vars; i++)
    {
        delete[] sep[i];
        delete[] sep2[i];
    }
    delete[] sep;
    delete[] sep2;
}

/////////////////////////////////////
/////////////////////////////////////


void pc_superset(int target, int *pc, int **sep, vector<int> &or_pc, vector<int> &and_pc)
{
    int i, j, last_added, stop, max_dep_index, remaining_length;
    int *sep2;
    double *dep, *cache_dep;
    double max_dep;
    MAP_DOUBLE_INT or_map, and_map;

    or_pc.clear();
    and_pc.clear();
    or_map.clear();
    and_map.clear();

    if (pc_superset_computed[target] == true)
    {
        for (i = 0; i < n_vars; i++)
        {
            pc[i] = 0;
        }
        for (vector<int>::iterator iter = pc_superset_cache[target].begin(); iter < pc_superset_cache[target].end(); iter++)
        {
            pc[*iter] = 1;
        }
        for (i = 0; i < n_vars; i++)
        {
            for (j = 0; j < max_max_cond_size; j++)
            {
                sep[i][j] = sep_cache[target][i][j];
            }
        }
        for (vector<int>::iterator iter = or_pc_cache[target].begin(); iter < or_pc_cache[target].end(); iter++)
        {
            or_pc.push_back(*iter);
        }
        for (vector<int>::iterator iter = and_pc_cache[target].begin(); iter < and_pc_cache[target].end(); iter++)
        {
            and_pc.push_back(*iter);
        }
        return;
    }

    for (i = 0; i < n_vars; i++)
        for (j = 0; j < max_max_cond_size; j++)
            sep[i][j] = -1;

    sep2 = new int[max_max_cond_size];

    dep = new double[n_vars];
    cache_dep = new double[n_vars];

    for (i = 0; i < n_vars; i++)
    {
        cache_dep[i] = (double)0.0;
    }

    for (i = 0; i < n_vars; i++)
        if (n_states[i] > 1 && n_states[target] > 1)
            pc[i] = 0;
        else
            pc[i] = -1;
    pc[target] = -1;

    last_added = -1;

    do
    {
        stop = 1;

        for (i = 0; i < n_vars; i++)
            dep[i] = (double)0.0;

        for (i = 0; i < n_vars; i++)
            if (pc[i] == 0)
            {
                dep[i] = min_dep_fast(i, target, last_added, pc, sep2);

                if (dep[i] <= (double)(-1.0) && last_added != -1)
                {
                    or_map[dep[i]] = i;
                }

                if (dep[i] == (double)0.0)
                    dep[i] = cache_dep[i];
                else
                {
                    if (cache_dep[i] == (double)0.0 || dep[i] < cache_dep[i])
                    {
                        cache_dep[i] = dep[i];

                        if (dep[i] <= (double)(-1.0))
                            for (j = 0; j < max_max_cond_size; j++)
                                sep[i][j] = sep2[j];
                    }
                    else
                        dep[i] = cache_dep[i];
                }
            }

        for (i = 0; i < n_vars; i++)
            if (dep[i] <= (double)(-1.0))
                pc[i] = -1;
            else if (dep[i] >= (double)1.0)
                stop = 0;

        if (stop == 0)
        {
            //last_added = k_greedy(dep);
            max_dep = 0.0;
            max_dep_index = -1;
            for (i = 0; i < n_vars; i++)
            {
                if (dep[i] > max_dep)
                {
                    max_dep = dep[i];
                    max_dep_index = i;
                }
            }
            last_added = max_dep_index;
            
            and_map[max_dep] = last_added;

            pc[last_added] = 1;
        }
        else
        {
            for (i = 0; i < n_vars; i++)
                dep[i] = (double)0.0;

            for (i = 0; i < n_vars; i++)
                if (pc[i] == 1 && i != last_added) //the last node added will not be cut.
                {
                    pc[i] = 0;

                    dep[i] = min_dep_fast(i, target, last_added, pc, sep2);

                    pc[i] = 1;

                    if (dep[i] <= (double)(-1.0))
                        for (j = 0; j < max_max_cond_size; j++)
                            sep[i][j] = sep2[j];
                }

            for (i = 0; i < n_vars; i++)
                if (pc[i] == 1 && i != last_added && dep[i] <= (double)(-1.0))
                {
                    pc[i] = -1;

                    for (MAP_DOUBLE_INT::iterator iter = and_map.begin(); iter != and_map.end(); iter++)
                    {
                        if (iter->second == i)
                        {
                            and_map.erase(iter);
                        }
                    }
                    
                    or_map[dep[i]] = i;

                    stop = 0;
                }
        }
    } while (stop == 0);

    for (i = 0; i < n_vars; i++)
        if (pc[i] == -1)
            pc[i] = 0;

    remaining_length = floor((and_map.size() * K_and) + 0.5);
    for (MAP_DOUBLE_INT::iterator iter = and_map.begin(); iter != and_map.end(); iter++)
    {
        if (remaining_length > 0)
        {
            and_pc.push_back(iter->second);
            and_pc_cache[target].push_back(iter->second);
            remaining_length--;
        }
    }

    remaining_length = floor((or_map.size() * K_or) + 0.5);
    for (MAP_DOUBLE_INT::reverse_iterator iter = or_map.rbegin(); iter != or_map.rend(); iter++)
    {
        if (remaining_length > 0)
        {
            or_pc.push_back(iter->second);
            or_pc_cache[target].push_back(iter->second);
            remaining_length--;
        }
    }

    pc_superset_cache_count++;
    pc_superset_computed[target] = true;
    for (i = 0; i < n_vars; i++)
        if (pc[i] == 1)
            pc_superset_cache[target].push_back(i);
    for (i = 0; i < n_vars; i++)
    {
        for (j = 0; j < max_max_cond_size; j++)
        {
            sep_cache[target][i][j] = sep[i][j];
        }
    }

    delete[] sep2;
    delete[] dep;
    delete[] cache_dep;
}

/////////////////////////////////////

double min_dep(int var, int target, int compulsory, int *pc, int *sep)
//It computes the minimum of the depence measure between var and target conditioned on
//subsets of pc[x1] st compulsory is always included (-1=no compulsory). Note that only
//those subsets of pc[x1] that have a size no larger than max_cond_size are considered.
//It returns the minimum dependence [-inf,-1] for H0 (the smaller the more independent),
//[1,+inf] for H1 (the higher the more dependent), or 0 for not enough data (see N_CASES_PER_DF),
//and sep[-1i].
{
    int i, j, n_pc, stop, cond_size, max_cond_size;
    double dep, aux;
    int *cond, *code, *cond_index;

    cond = new int[max_max_cond_size];
    for (i = 0; i < max_max_cond_size; i++)
    {
        cond[i] = -1;
        sep[i] = -1;
    }

    n_pc = 0;
    for (i = 0; i < n_vars; i++)
        if (pc[i] == 1 && i != compulsory)
            n_pc++;

    max_cond_size = max_max_cond_size;
    if (compulsory != -1)
        max_cond_size--;
    if (max_cond_size > n_pc)
        max_cond_size = n_pc;

    if (k_conditon != -1)
    {
        if (compulsory != -1)
        {
            if (max_cond_size >= k_conditon)
            {
                max_cond_size = k_conditon - 1;
            }
        }
        else
        {
            if (max_cond_size > k_conditon)
            {
                max_cond_size = k_conditon;
            }
        }
    }

    code = new int[n_pc];
    j = 0;
    for (i = 0; i < n_vars; i++)
        if (pc[i] == 1 && i != compulsory)
        {
            code[j] = i;

            j++;
        }

    dep = (double)0.0;
    for (cond_size = 0; cond_size <= max_cond_size; cond_size++)
    {
        cond_index = new int[cond_size];
        for (i = 0; i < cond_size; i++)
            cond_index[i] = i;

        do
        {
            stop = 0;

            for (i = 0; i < cond_size; i++)
                cond[i] = code[cond_index[i]];
            if (compulsory != -1)
                cond[cond_size] = compulsory;

            aux = compute_dep(var, target, cond);

            if (aux != (double)0.0 && (dep == (double)0.0 || aux < dep))
            {
                dep = aux;

                if (dep <= (double)(-1.0))
                {
                    for (i = 0; i < max_max_cond_size; i++)
                        sep[i] = cond[i];

                    stop = 1;
                    cond_size = max_cond_size + 1; //end loop early
                }
            }

            if (stop == 0)
                stop = next_cond_index(n_pc, cond_size, cond_index);
        } while (stop == 0);

        delete[] cond_index;
    }

    delete[] code;
    delete[] cond;

    return dep;
}

/////////////////////////////////////

double min_dep_fast(int var, int target, int compulsory, int *pc, int *sep)
//It computes the minimum of the depence measure between var and target conditioned on
//subsets of pc[x1] st compulsory is always included (-1=no compulsory). Note that only
//those subsets of pc[x1] that have a size no larger than max_cond_size are considered.
//It returns the minimum dependence [-inf,-1] for H0 (the smaller the more independent),
//[1,+inf] for H1 (the higher the more dependent), or 0 for not enough data (see N_CASES_PER_DF),
//and sep[-1i].
{
    int i, j, n_pc, stop, cond_size, max_cond_size;
    double dep, aux;
    int *cond, *code, *cond_index;

    cond = new int[max_max_cond_size];
    for (i = 0; i < max_max_cond_size; i++)
    {
        cond[i] = -1;
        sep[i] = -1;
    }

    n_pc = 0;
    for (i = 0; i < n_vars; i++)
        if (pc[i] == 1 && i != compulsory)
            n_pc++;

    max_cond_size = max_max_cond_size;
    if (compulsory != -1)
        max_cond_size--;
    if (max_cond_size > n_pc)
        max_cond_size = n_pc;

    if (k_conditon != -1)
    {
        if (compulsory != -1)
        {
            if (max_cond_size >= k_conditon)
            {
                max_cond_size = k_conditon - 1;
            }
        }
        else
        {
            if (max_cond_size > k_conditon)
            {
                max_cond_size = k_conditon;
            }
        }
    }

    code = new int[n_pc];
    j = 0;
    for (i = 0; i < n_vars; i++)
        if (pc[i] == 1 && i != compulsory)
        {
            code[j] = i;

            j++;
        }

    dep = (double)0.0;
    for (cond_size = 0; cond_size <= max_cond_size; cond_size++)
    {
        cond_index = new int[cond_size];
        for (i = 0; i < cond_size; i++)
            cond_index[i] = i;

        do
        {
            stop = 0;

            for (i = 0; i < cond_size; i++)
                cond[i] = code[cond_index[i]];
            if (compulsory != -1)
                cond[cond_size] = compulsory;

            aux = Search_dep(var, target, cond);

            if (aux != (double)0.0 && (dep == (double)0.0 || aux < dep))
            {
                dep = aux;

                if (dep <= (double)(-1.0))
                {
                    for (i = 0; i < max_max_cond_size; i++)
                        sep[i] = cond[i];

                    stop = 1;
                    cond_size = max_cond_size + 1; //end loop early
                }
            }

            if (stop == 0)
                stop = next_cond_index(n_pc, cond_size, cond_index);
        } while (stop == 0);

        delete[] cond_index;
    }

    delete[] code;
    delete[] cond;

    return dep;
}

/////////////////////////////////////

int next_cond_index(int n_pc, int cond_size, int *cond_index)
//Non-recursive subset function
{
    int i, j, stop;

    stop = 1;
    for (i = cond_size - 1; i >= 0; i--)
    {
        if (cond_index[i] < n_pc + i - cond_size)
        {
            cond_index[i]++;

            if (i < cond_size - 1)
                for (j = i + 1; j < cond_size; j++)
                    cond_index[j] = cond_index[j - 1] + 1;

            stop = 0;
            i = -1;
        }
    }

    return stop;
}

/////////////////////////////////////

double compute_dep(int var, int target, int *cond)
//It computes the depence measure between var and target conditioned on cond[-1i], where
//cond[-1i] is an array of size max_max_cond_size.
//It returns the dependence [-inf,-1] for H0 (the smaller the more independent), [+1,+inf]
//for H1 (the higher the more dependent), or 0 for not enough data (see N_CASES_PER_DF).

//It can be speeded up by allocating memory only once in main() at the cost of memory space,
//and by caching the critical value of the test as a function of the degrees of freedom.
//Only ss[][][] is strictly necessary. The rest of the arrays are used to gain speed at the
//cost of memory space.
{
    int i, j, k, n_cond_states, cond_state, df, df_target, df_var;
    double statistic, dep, aux;
    int *ss_cond;
    int **ss_target, **ss_var;
    int ***ss;

    average_n_G2++;

    n_cond_states = 1;
    for (i = 0; i < max_max_cond_size; i++)
        if (cond[i] != -1)
            n_cond_states = n_cond_states * n_states[cond[i]];

    if (n_cond_states * (n_states[target] - 1) * (n_states[var] - 1) * (int)N_CASES_PER_DF > n_cases)
        return (double)0.0;

    ss_cond = new int[n_cond_states];
    for (i = 0; i < n_cond_states; i++)
        ss_cond[i] = 0;

    ss_target = new int *[n_cond_states];
    for (i = 0; i < n_cond_states; i++)
    {
        ss_target[i] = new int[n_states[target]];
        for (j = 0; j < n_states[target]; j++)
            ss_target[i][j] = 0;
    }

    ss_var = new int *[n_cond_states];
    for (i = 0; i < n_cond_states; i++)
    {
        ss_var[i] = new int[n_states[var]];
        for (j = 0; j < n_states[var]; j++)
            ss_var[i][j] = 0;
    }

    ss = new int **[n_cond_states];
    for (i = 0; i < n_cond_states; i++)
    {
        ss[i] = new int *[n_states[target]];
        for (j = 0; j < n_states[target]; j++)
        {
            ss[i][j] = new int[n_states[var]];
            for (k = 0; k < n_states[var]; k++)
                ss[i][j][k] = 0;
        }
    }

    for (i = 0; i < n_cases; i++)
    {
        cond_state = 0;
        for (j = 0; j < max_max_cond_size; j++)
            if (cond[j] != -1)
                cond_state = cond_state * n_states[cond[j]] + data_processed[i][cond[j]];

        ss_cond[cond_state]++;
        ss_target[cond_state][data_processed[i][target]]++;
        ss_var[cond_state][data_processed[i][var]]++;
        ss[cond_state][data_processed[i][target]][data_processed[i][var]]++;
    }

    //X^2 statistic.
    /*
    statistic=(double)0.0;
    for(i=0;i<n_cond_states;i++)
    if(ss_cond[i]>0)
    for(j=0;j<n_states[target];j++)
    if(ss_target[i][j]>0)
    for(k=0;k<n_states[var];k++)
    if(ss_var[i][k]>0)
    {aux=(double)(ss_target[i][j]*ss_var[i][k])/(double)ss_cond[i];
    statistic=statistic+pow((double)ss[i][j][k]-aux,(double)2.0)/aux;
    }
    */

    //G^2 statistic based on mutual information.
    /*
    statistic=(double)0.0;
    for(i=0;i<n_cond_states;i++)
    if(ss_cond[i]>0)
    for(j=0;j<n_states[target];j++)
    if(ss_target[i][j]>0)
    for(k=0;k<n_states[var];k++)
    if(ss[i][j][k]>0)
    statistic=statistic+(double)ss[i][j][k]*(log10((double)ss[i][j][k])-
    log10((double)ss_target[i][j])-log10((double)ss_var[i][k])+
    log10((double)ss_cond[i]))/log10((double)2.0);
    statistic=statistic*(double)2.0;
    */

    //G^2 statistic based on observed and expected frequencies.

    statistic = (double)0.0;
    for (i = 0; i < n_cond_states; i++)
        if (ss_cond[i] > 0)
            for (j = 0; j < n_states[target]; j++)
                if (ss_target[i][j] > 0)
                    for (k = 0; k < n_states[var]; k++)
                        if (ss[i][j][k] > 0)
                        {
                            aux = (double)(ss_target[i][j] * ss_var[i][k]) / (double)ss_cond[i];
                            statistic = statistic + (double)ss[i][j][k] * (log((double)ss[i][j][k]) - log(aux));
                        }
    statistic = statistic * (double)2.0;

    //Reduced df due to zero entries. This does not make sense with the X^2 statistic.
    /*
    df=n_cond_states*(n_states[target]-1)*(n_states[var]-1);
    for(i=0;i<n_cond_states;i++)
    for(j=0;j<n_states[target];j++)
    for(k=0;k<n_states[var];k++)
    if(ss[i][j][k]==0)
    df--;
    */

    //Reduced df due to zero marginals.

    df = 0;
    for (i = 0; i < n_cond_states; i++)
        if (ss_cond[i] > 0)
        {
            df_target = 0;
            for (j = 0; j < n_states[target]; j++)
                if (ss_target[i][j] > 0)
                    df_target++;

            df_var = 0;
            for (k = 0; k < n_states[var]; k++)
                if (ss_var[i][k] > 0)
                    df_var++;

            df = df + (df_target - 1) * (df_var - 1);
        }

    delete[] ss_cond;
    for (i = 0; i < n_cond_states; i++)
    {
        delete[] ss_target[i];
        delete[] ss_var[i];
        for (j = 0; j < n_states[target]; j++)
            delete[] ss[i][j];
        delete[] ss[i];
    }
    delete[] ss_target;
    delete[] ss_var;
    delete[] ss;

    if (statistic < (double)0.0) //statistic sometimes takes values -0.0. Due to loss of precision ?
        statistic = (double)0.0;

    if (df <= 0) //Naive ?
        df = 1;

    dep = gammq((double)0.5 * (double)df, (double)0.5 * statistic);

    if (dep < (double)0.0) //Just in case.
        dep = (double)0.0;
    else if (dep > (double)1.0)
        dep = (double)1.0;

    if (dep <= alpha)
    {
        dep = (double)2.0 - dep;

        if (dep == (double)2.0)
            dep = (double)2.0 + statistic / (double)df;

        return dep;
    }
    else
    {
        dep = (double)(-1.0) - dep;

        if (dep == (double)(-2.0))
            dep = (double)(-2.0) - statistic / (double)df;

        return dep;
    }
}

/////////////////////////////////////

double Search_dep(int var, int target, int *cond)
{
    int i, n_conds;
    std::string key_string, s_aux;
    std::stringstream ss;
    double dep;

    Search_Count2++;
    n_conds = 0;
    key_string = "";
    if (var < target)
    {
        ss << var;
        ss >> s_aux;
        key_string += s_aux;
        key_string += "_";
        ss.clear();
        ss << target;
        ss >> s_aux;
        key_string += s_aux;
    }
    else
    {
        ss << target;
        ss >> s_aux;
        key_string += s_aux;
        key_string += "_";
        ss.clear();
        ss << var;
        ss >> s_aux;
        key_string += s_aux;
    }

    for (i = 0; i < max_max_cond_size; i++)
    {
        if (cond[i] != -1)
        {
            n_conds++;
            key_string += "_";
            ss.clear();
            ss << cond[i];
            ss >> s_aux;
            key_string += s_aux;
        }
    }
    MAP_STRING_DOUBLE::iterator iter = Cond_Size[n_conds].find(key_string);
    if (iter != Cond_Size[n_conds].end())
    {
        Search_Count++;
        return iter->second;
    }

    dep = compute_dep(var, target, cond);
    Cond_Size[n_conds].insert(MAP_STRING_DOUBLE::value_type(key_string, dep));
    return dep;
}

/////////////////////////////////////

int compute_rows(const char *path)
{
    std::ifstream fileStream;
    std::string tmp;
    int count = 0;                       // 行数计数器
    fileStream.open(path, std::ios::in); //ios::in 表示以只读的方式读取文件
    if (fileStream.fail())               //文件打开失败:返回0
    {
        std::cout << "样本数据路径error" << std::endl;
        return 0;
    }
    else //文件存在
    {
        while (getline(fileStream, tmp, '\n')) //读取一行
        {
            if (tmp.size() > 0)
                count++;
        }
    }
    fileStream.close();
    return count;
}

/////////////////////////////////////

int compute_columns(const char *path)
{
    std::ifstream fileStream;
    fileStream.open(path, std::ios::in);
    if (!fileStream)
    {
        std::cout << "样本数据路径error" << std::endl;
        exit(0);
    }

    double tmp = 0;
    int count = 0; // 列数计数器
    char c;        //当前位置的字符
    c = fileStream.peek();
    while ((('\n' != c) && ('\r' != c)) && (!fileStream.eof())) // 指针指向的当前字符，仅观测，不移动指针位置
    {
        fileStream >> tmp;
        ++count;
        c = fileStream.peek();
    }

    fileStream.close();
    return count;
}

/////////////////////////////////////

void report_mb(int target, int *mb)
{
    std::fstream file;
    int i;

    average_perform_time = average_perform_time + ((double)(clock() - clock_init) / CLOCKS_PER_SEC);

    std::cout << "***Algorithm execution completed!***" << std::endl;
    std::cout << "Output:" << std::endl;

    //write to text
    file.open(mb_write_in_file_path, std::ios::out | std::ios::app); //read from memory into disk and open file in append mode
    if (!file)
    {
        std::cout << "error" << std::endl;
        return;
    }

    for (i = 0; i < n_vars; i++)
    {
        if (mb[i] == 1)
        {
            file << i+1 << " ";
        }
    }
    file << std::endl;
    file.close();


    std::cout << "mb (causal features) of label: ";
    for (i = 0; i < n_vars; i++)
    {
        if (mb[i] == 1)
        {
            std::cout << i+1 << " ";
        }
    }
    std::cout << std::endl;
}

/////////////////////////////////////