/*************************************************************************
 
 Travelling Salesman Problem with Time Windows
 
 ---------------------------------------------------------------------
 
 Copyright (c) 2009
 Manuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be>
 Christian Blum <cblum@lsi.upc.edu>
 
 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at: http://www.gnu.org/licenses/gpl.html
 
 ---------------------------------------------------------------------
 
 Compilation:
 
 g++ -o check_solution check_solution.cpp
 
 ---------------------------------------------------------------------
 
 References:
 
 [1] Manuel Lopez-Ibanez and Christian Blum. Beam-ACO for the
 travelling salesman problem with time windows. Computers &
 Operations Research, 37(9):1570–1583, 2010.
 doi:10.1016/j.cor.2009.11.015
 
 *************************************************************************/

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <cmath>
#include <cassert>
#include <string>
#include <algorithm>
#include <time.h>
using namespace std;

#define program_invocation_short_name "check_solution"

typedef double number_t;
#define STRING_NUMBER_IS "double"
#define NUMBER_T_MAX DBL_MAX
#define NUMBER_T_MIN DBL_MIN


class Solution {
    
public:
    
    static string instance;
    static int n;   // number of customers
    static void LoadInstance (string filename);
    
    vector<int> permutation;
    
    int _constraint_violations;
    
    Solution (string filename);
    
    void print_one_line (FILE *stream=stdout) const;
    
    number_t makespan() const;
    void evaluate();
    
    // private:
    
    number_t _makespan;
    number_t _tourcost; // Sum of the traversal cost along the tour.
    
    // time-window start
    static vector<number_t> window_start;
    
    // time-window end
    static vector<number_t> window_end;
    
    // travel time/distance
    static vector<vector<number_t> > distance;
    
};

string Solution::instance;

// number of customers
int Solution::n = 0;

// time-window start
vector<number_t> Solution::window_start;

// time-window end
vector<number_t> Solution::window_end;

// travel time/distance
vector<vector<number_t> > Solution::distance;


Solution::Solution(string filename)
: permutation (1,0), // Start at the depot.
_constraint_violations (0),
_makespan (0),
_tourcost (0)
{
    permutation.reserve (n+1);
    
    ifstream indata;
    int node;
    
    indata.open (filename.c_str());
    if (!indata) { // file couldn't be opened
        cout << "error: file " << filename.c_str() << " could not be opened"
        << endl;
        exit (EXIT_FAILURE);
    }
    
    for (int i = 1; i < n; i++) {
        indata >> node;
        permutation.push_back (node);
    }
    
    permutation.push_back (0); // Finish at the depot;
    
    if (!indata) {
        cout << "error: file " << filename.c_str()
        << " contains invalid solution" << endl;
        exit (EXIT_FAILURE);
    }
    
    indata.close();
}

void
Solution::LoadInstance (string filename)
{
    ifstream indata;
    number_t rtime;
    number_t ddate;
    
    instance = filename;
    indata.open (instance.c_str());
    if (!indata) { // file couldn't be opened
        cout << "error: file " << instance.c_str() << " could not be opened"
        << endl;
        exit (EXIT_FAILURE);
    }
    
    // Customer 0 is the depot.
    indata >> n;
    if (!indata || n <= 0) {
        cout << "error: invalid number of customers" << endl;
        exit (EXIT_FAILURE);
    }
    
    window_start.reserve(n);
    window_end.reserve(n);
    distance.reserve(n);
    
    for (int i = 0 ; i < n; i++) {
        distance.push_back(vector<number_t>(n,0));
        for (int j = 0; j < n; j++) {
            indata >> distance[i][j];
        }
    }
    
    if (!indata) {
        cout << "error: invalid time windows" << endl;
        exit (EXIT_FAILURE);
    }
    
    for (int i = 0; i < n; i++) {
        indata >> rtime >> ddate;
        window_start.push_back (rtime);
        window_end.push_back(ddate);
    }
    
    if (!indata) {
        cout << "error: invalid distance matrix" << endl;
        exit (EXIT_FAILURE);
    }
    
    indata.close();
}

void Solution::print_one_line (FILE *stream) const
{
    fprintf (stream, "makespan = %.2f\ttourcost = %.2f\tconstraint violations = %d\tpermutation =",
             double (makespan()), double(_tourcost), _constraint_violations);
    
    // Customers 0 and n+1 are always the depot (permutation[0] == 0,
    // permutation[n] == 0), so do not print them.
    for (int i = 1; i < int(permutation.size()) - 1; i++) {
        fprintf (stream, " %d", permutation[i]);
    }
    fprintf (stream, "\n");
}

number_t Solution::makespan() const
{
    return _makespan;
}

void Solution::evaluate()
{
    _makespan = 0;
    _tourcost = 0;
    int prev = 0; // starts at the depot
    int cviols = 0;
    
    if (int(permutation.size() - 1) != n) {
        printf ("invalid: (permutation.size() == %d) != (n == %d)\n",
                int(permutation.size() - 1), n);
        exit (EXIT_FAILURE);
    }
    
    for (int i = 1; i < n; i++) {
        int node = permutation[i];
        
        _tourcost += distance[prev][node];
        _makespan = max (_makespan + distance[prev][node], window_start[node]);
        
        if (_makespan > window_end[node]) {
            cviols++;
        }
        prev = node;
    }
    
    // finish at the depot
    _tourcost += distance[prev][0];
    
    _makespan = max (_makespan + distance[prev][0], window_start[0]);
    
    if (_makespan > window_end[0])
        cviols++;
    
    _constraint_violations = cviols;
}

static void version(void)
{
    printf ("%s", program_invocation_short_name);
    printf("\n\n"
           "Copyright (C) 2009\n"
           "Manuel Lopez-Ibanez (manuel.lopez-ibanez@ulb.ac.be) and\n"
           "Christian Blum (cblum@lsi.upc.edu)\n"
           "\n"
           "This is free software, and you are welcome to redistribute it under certain\n"
           "conditions.  See the GNU General Public License for details. There is NO   \n"
           "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"
           "\n"        );
}

static void usage(void)
{
    printf("\n"
           "Usage: %s INSTANCE_FILE SOLUTION_FILE\n\n", program_invocation_short_name);
    
    printf(
           "Reads an instance file and a file with a permutation (from 1 to N, that is, not containing the depot) and evaluates the solution.\n"
           "\n");
    version ();
}

// algorithm of using VNS_1_OPT with forward and backward
// and destruct and construct way of heuristics
Solution VIG_VNS(Solution s0);
Solution VNS_1_OPT(Solution s);
Solution DestructConstruct(Solution s, int d);
bool lex_s1_greater_than_s2(Solution s1, Solution s2);
Solution Backward_1_OPT(Solution s);
Solution Forward_1_OPT(Solution s);
bool isFeasible(int i, int j, Solution s);
Solution insert(Solution s, int i, int j);

int main(int argc, char **argv)
{
    if (argc != 3) {
        printf ("error: invalid number of parameters\n");
        usage();
        exit (EXIT_FAILURE);
    }
    string input_instance = argv[1];
    string input_solution = argv[2];
    
    Solution::LoadInstance (input_instance);
    Solution s(input_solution);
    s.evaluate();
    s.print_one_line();
    
    s = VIG_VNS(s);
    
    s.evaluate();
    s.print_one_line();
    printf ("\n");
    
    return EXIT_SUCCESS;
}

Solution VIG_VNS(Solution s0)
{
    int n = s0.n - 1;
    int kmax = floor((n - 1)/5);
    Solution s = VNS_1_OPT(s0);
    Solution s_best = s;
    s_best.evaluate();
    cout << "does VNS_1_OPT really work" << " ";
    s_best.print_one_line();
    
    // decide to use number of iterations
    int index = 0;
    while (index < 80000) {
        // cout << "------------In loop " << index << "--------------" << endl;
        int k = 1;
        
        while (k <= kmax) {
            // cout << "k = " << k << endl;
            int d = k*5;
            Solution s1 = DestructConstruct(s, d);
            Solution s2 = VNS_1_OPT(s1);
            
            s1.evaluate();
            s2.evaluate();
            // cout << "After construct and opt " << endl;
            // s1.print_one_line();
            // s2.print_one_line();
            
            //            if (lex_s1_greater_than_s2(s, s2))
            if (lex_s1_greater_than_s2(s2, s))
            {
                k = 1;
                s = s2;
                //                if (lex_s1_greater_than_s2(s_best, s))
                if (lex_s1_greater_than_s2(s, s_best))
                {
                    s_best = s;
                    s_best.evaluate();
                    s_best.print_one_line();
                }
            }
            else {
                k = k + 1;
            }
        }
        
        // s_best.evaluate();
        // cout << "At the end of loop " << index << " ";
        // s_best.print_one_line();
        
        index++;
    }
    
    return s_best;
}

Solution VNS_1_OPT(Solution s)
{
    Solution choice = s;
    int kmax = 2;
    int k = 1;
    while(k <= kmax) {
        Solution s1 = s;
        
        if (k == 1)
            s1 = Backward_1_OPT(s);
        if (k == 2)
            s1 = Forward_1_OPT(s);
        
        s1.evaluate();
        // cout << "s1.constraints_voi = " << s1._constraint_violations << " s.constraints_voi = " << s._constraint_violations << endl;
        if(lex_s1_greater_than_s2(s1, s)) {
            //        if(lex_s1_greater_than_s2(s, s1)) {
            k = 1;
            // choice = s1;
            s = s1;
            // cout << "--------**************#######" << endl;
        }
        else {
            k = k + 1;
        }
    }
    return s;
}

Solution Backward_1_OPT(Solution s)
{
    Solution old_solution = s;
    // for(int i = s.n - 2; i >= 2; i--) {
    for(int i = s.n - 1; i >= 2; i--) {
       int remove = s.distance[i-1][i] + s.distance[i][i+1] + s.distance[i-1][i+1];
        // int remove = s.distance[s.permutation[i-1]][s.permutation[i]] +
        //              s.distance[s.permutation[i]][s.permutation[i+1]] +
        //              s.distance[s.permutation[i-1]][s.permutation[i+1]];
        // int remove = s.distance.at(i-1).at(i);
        // int remove = s.distance[i-1][i];
        
        // erase the ith element (1,2,3,4, .....)
        // int node_removed = s.permutation.at(i-1);
        //        s.permutation.erase(s.permutation.begin()+i-1);
        //        cout << "new slotion size = " << s.permutation.size() << endl;
        //        cout << "old slotion size = " << old_solution.permutation.size() << endl;
        // s.print_one_line();
        
        for (int j = i-1; j >= 1; j--)
        {
           int add = s.distance[j-1][i] + s.distance[i][j] + s.distance[j-1][j];
            // int add = s.distance[s.permutation[j-1]][s.permutation[i]] +
            // s.distance[s.permutation[i]][s.permutation[j]] +
            // s.distance[s.permutation[j-1]][s.permutation[j]];
            int gain = add - remove;
            
            // cout << "gain = " << gain << endl;
            if (gain < 0 && isFeasible(i, j, s)) // and is feasible (i, j) ???????????????
            {
                //                s.print_one_line();
                //                s.permutation.insert(s.permutation.begin() + j - 1, node_removed);
                //                cout << "new slotion size = " << s.permutation.size() << endl;
                //                s.print_one_line();
                Solution s1 = insert(s, i, j);
                
                //                if (lex_s1_greater_than_s2(old_solution, s))
                if (lex_s1_greater_than_s2(s1, s))
                {
                    // return s;
                    // cout << "Backword works~~~~~~~~~~~~~~~~" << endl;
                    s = s1;
                }
            }
        }
        // s = old_solution;
    }
    return s;
}

Solution insert(Solution s, int i, int j)
{
    Solution s1 = s;
    int node_removed = s1.permutation.at(i);
    s1.permutation.erase(s1.permutation.begin()+i);
    s1.permutation.insert(s1.permutation.begin() + j, node_removed);
    return s1;
}

Solution Forward_1_OPT(Solution s)

{
    Solution old_solution = s;
    // for(int i = 2; i <= s.n - 2; i++) {
    for(int i = 1; i <= s.n - 2; i++) {
       int remove = s.distance[i-1][i] + s.distance[i][i+1] + s.distance[i-1][i+1];
        // int remove = s.distance[s.permutation[i-1]][s.permutation[i]] +
        // s.distance[s.permutation[i]][s.permutation[i+1]] +
        // s.distance[s.permutation[i-1]][s.permutation[i+1]];
        
        // erase the ith element (1,2,3,4, .....)
        int node_removed = s.permutation.at(i-1);
        //        s.permutation.erase(s.permutation.begin()+i-1);
        
        for (int j = i+1; j <= s.n - 1; j++)
        {
           int add = s.distance[i][j+1] + s.distance[j][i] + s.distance[j][j+1];
            // int add = s.distance[s.permutation[i]][s.permutation[j+1]] +
            // s.distance[s.permutation[j]][s.permutation[i]] +
            // s.distance[s.permutation[j]][s.permutation[j+1]];
            int gain = add - remove;
            if (gain < 0 && isFeasible(i, j, s)) // and is feasible (i, j) ???????????????
            {
                //                s.permutation.insert(s.permutation.begin() + j - 1, node_removed);
                Solution s1 = insert(s, i, j);
                //                if (lex_s1_greater_than_s2(old_solution, s))
                if (lex_s1_greater_than_s2(s1, s))
                {
                    // return s;
                    // cout << "Forword works~~~~~~~~~~~~~~~~" << endl;
                    s = s1;
                }
            }
        }
    }
    return s;
}

// lexicographic order s1 is greater than s2
// 1. s1 feasible, s2 infeasible
// 2. s1 has smaller objective, if both feasible
// 3. s1 has less number of violations, if both infeasbile
bool lex_s1_greater_than_s2(Solution s1, Solution s2) {
    s1.evaluate();
    s2.evaluate();
    
    if (s1._constraint_violations == 0 && s2._constraint_violations != 0)
        return true;
    
    if (s1._constraint_violations == 0 && s2._constraint_violations == 0) {
        if (s1._tourcost < s2._tourcost)
            return true;
    }
    
    if (s1._constraint_violations != 0 && s2._constraint_violations != 0) {
        if (s1._constraint_violations < s2._constraint_violations)
            return true;
    }
    
    return false;
}

bool isFeasible(int i, int j, Solution s) {
    if (s.window_start[j] + s.distance[j][i] <= s.window_end[i] &&
        s.window_start[i] + s.distance[i][j+1] <= s.window_end[j+1])
        return true;
    return false;
}

// int transformIfInB(const std::vector<int>& B, int ptr)
// {
//   return find(B.begin(),B.end(), ptr) != B.end() ? 0 : ptr;
// }
// void shuffle(int *arr, size_t n)
// {
//     if (n > 1)
//     {
//         size_t i;
//         srand(time(NULL));
//         for (i = 0; i < n - 1; i++)
//         {
//           size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
//           int t = arr[j];
//           arr[j] = arr[i];
//           arr[i] = t;
//         }
//     }
// }


Solution DestructConstruct(Solution s, int d) {
    Solution old_solution = s;
    int max = s.n - 2;
    int min = 1;
    std::vector<int> removed;
    std::vector<int> deducted;
    std::vector<int> origin = s.permutation;
    
    removed.clear();
    
    //    cout << "d = " << d << endl;
    
    // get vector of removed elements randomly choose
    time_t t;
    
    srand(time(&t));
    // srand(time(0)); // Seed the time
    while(removed.size() <= d) {
        // Generate the number, assign to variable.
        int index = rand()%(max-min)+min;
        
        int target = old_solution.permutation[index];
        bool isPresent = (std::find(removed.begin(), removed.end(), target) != removed.end());
        if (!isPresent)
            removed.push_back(target);
    }
    
    //    for (int i = 0; i < removed.size(); ++i)
    //    {
    //        cout << removed[i] << " ";
    //    }
    //    cout << endl;
    
    // get solution without removed elements
    for (int i = 0; i < removed.size(); ++i)
    {
        origin.erase(std::remove(origin.begin(), origin.end(), removed[i]), origin.end());
    }
    
    
    // cout << "origin = " << endl;
    //    for (int i = 0; i < origin.size(); ++i)
    //    {
    //        cout << origin[i] << " ";
    //    }
    //    cout << endl;
    
    for (int i = 0; i < removed.size(); ++i)
    {
        // add one element back to orgin
        int position = origin.size() + 1;
        int best_postion_cost = 0;
        int best_postion = 0;
        
        // find best postion
        for (int j = 1; j < position-1; ++j)
        {
            std::vector<int> tmp = origin;
            
            int tour_cost = 0;
            
            tmp.insert(tmp.begin() + j, removed[i]);
            
            int m = 0;
            while (m < origin.size())
            {
                tour_cost = tour_cost + s.distance[tmp[m]][tmp[m+1]];
                m++;
            }
            
            if(j == 1) {
                best_postion_cost = tour_cost;
                best_postion = 1;
            }
            
            else if(tour_cost < best_postion_cost) {
                best_postion_cost = tour_cost;
                best_postion = j;
            }
        }
        //        cout << "best postion is " << best_postion<< " insert " << removed[i] << endl;
        origin.insert(origin.begin() + best_postion, removed[i]);
        
        //        for (int i = 0; i < origin.size(); ++i)
        //        {
        //            cout << origin[i] << " ";
        //        }
        //        cout << endl;
    }
    
    //    cout << "best postion finded:" << endl;
    //    for (int i = 0; i < origin.size(); ++i)
    //    {
    //        cout << origin[i] << " ";
    //    }
    s.permutation = origin;
    s.evaluate();
    
    int NFT_new = s._constraint_violations^10;
    int NFT_old = old_solution._constraint_violations^10;
    int obj_new = s._tourcost + NFT_new;
    int obj_old = old_solution._tourcost + NFT_old;
    
    return s;
    // if (obj_new < obj_old)
    //     return s;
    // else return old_solution;
}
