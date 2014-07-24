/*
****************************** Copyright & License *****************************
CALIPER v1.0 is a software framework for the reliability lifetime evaluation of Multicore architectures. Copyright (C) 2014 Politecnico di Milano.

This framework is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This framework is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details (http://www.gnu.org/licenses/).

Neither the name of Politecnico di Milano nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
********************************************************************************
*/

#include <iostream>
#include <cstdlib>
#include <map>
#include "utils.h"

/*
 * Tokenize a string on the basis of a delimiter
 */
void stringTokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {
    tokens.clear();
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

/*
 * Load operating configurations from a file
 */
void loadConfigurations(const std::string& fileName, std::map<std::string, std::vector<std::vector<double> > >& configurations, int& init_num_of_cores,
        int& min_num_of_cores) {
    //load configurations info from file
    int lineNumber = 0, i;
    std::string line;
    std::vector<std::string> tokens, deads;
    std::string config;
    std::vector<double> dagings;
    bool initNotRead = true;

    std::ifstream confFile(fileName.c_str());
    if (confFile.is_open()) {
        getline(confFile, line);
        lineNumber++;
        while (!confFile.eof()) {
            if (line != "") {
                //read current configuration
                stringTokenize(line, tokens, " ");
                config = tokens[0];
                if (tokens.size() < 2) {
                    std::cerr << "Line " << lineNumber << " - Missing data" << std::endl;
                    exit(1);
                }
                for (i = 1; i < tokens.size(); i++) {
                    dagings.push_back(atof(tokens[i].c_str()));
                }
                //REMOVED TO SUPPORT MULTIPLE PERIODIC CONFIGURATIONS
                //if (configurations.count(config) != 0) {
                //  std::cerr << "Line " << lineNumber << " - Configuration already specified: " << config << std::endl;
                //  exit(1);
                //}
                configurations[config].push_back(dagings);
                dagings.clear();

                //update min and max number of cores
                if (config == EMPTY_SET) {
                    init_num_of_cores = tokens.size() - 1; //-1 since we have to discard the configuration string
                    initNotRead = false;
                    min_num_of_cores = init_num_of_cores;
                } else {
                    if (initNotRead) {
                        std::cerr << "Line " << lineNumber << " - The first line of the file must be the initial situation" << std::endl;
                        exit(1);
                    }
                    stringTokenize(config, deads, ",");
                    //check of the current configuration
                    std::set<int> checksum;
                    for (int j = 0; j < deads.size(); j++) {
                        int currCore = atoi(deads[j].c_str());
                        if (currCore >= init_num_of_cores || currCore < 0) {
                            std::cerr << "Line " << lineNumber << " - Invalid core number: " << config << std::endl;
                            exit(1);
                        }
                        if (checksum.find(currCore) != checksum.end()) {
                            std::cerr << "Line " << lineNumber << " - Not valid sequence of deads: " << config << std::endl;
                            exit(1);
                        }
                        checksum.insert(currCore);
                    }
                    if (init_num_of_cores - deads.size() < min_num_of_cores)
                        min_num_of_cores = init_num_of_cores - deads.size();
                }
            }
            getline(confFile, line);
            lineNumber++;
        }
        confFile.close();
    } else {
        std::cerr << "Configuration file not found: " << fileName << std::endl;
        exit(1);
    }
}

/*
 * Functions for the generation of the tree
 */
void reverse(int *ar, int len);
int MBnext_permutation(int *ar, int len);
int MBnext_k_permutation(int *ar, int n, int k);
int MBnext_combination(int *ar, int n, int k);

void generateTree(int max, int min, std::set<std::string>& generated) {
    int *numbers;
    int i, k, j;
    std::string currConf;

    //generate root
    generated.insert(EMPTY_SET);

    //generate first level
    for (i = 0; i < max; i++)
        generated.insert(MAKE_STRING(i));

    //generate other levels
    numbers = new int[max - min + 1];
    for (k = 2; k <= max - min; k++) {

        for (i = 0; i < k; i++)
            numbers[i] = i;
        do {
            currConf = MAKE_STRING(numbers[0]);
            for (j = 1; j < k; j++)
                currConf = currConf + "," + MAKE_STRING(numbers[j]);
            generated.insert(currConf);
        } while (MBnext_k_permutation(numbers, max, k));
    }
    delete[] numbers;
}

void reverse(int *ar, int len) {
    int i, j, tmp;
    for (i = 0, j = len - 1; i < j; i++, j--) {
        tmp = ar[i];
        ar[i] = ar[j];
        ar[j] = tmp;
    }
}

int MBnext_permutation(int *ar, int len) {
    int i1, i2, tmp;
    int result = 0;

    /* Find the rightmost element that is the first in a pair in ascending order */
    for (i1 = len - 2, i2 = len - 1; ar[i2] <= ar[i1] && i1 != 0; i1--, i2--)
        ;
    if (ar[i2] <= ar[i1]) {
        /* If not found, array is highest permutation */
        reverse(ar, len);
    } else {
        /* Find the rightmost element to the right of i1 that is greater than ar[i1] */
        for (i2 = len - 1; i2 > i1 && ar[i2] <= ar[i1]; i2--)
            ;
        /* Swap it with the first one */
        tmp = ar[i1];
        ar[i1] = ar[i2];
        ar[i2] = tmp;

        /* Reverse the remainder */
        reverse(ar + i1 + 1, len - i1 - 1);
        result = 1;
    }
    return result;
}

int MBnext_k_permutation(int *ar, int n, int k) {
    int result = MBnext_permutation(ar, k);
    if (result == 0) {
        result = MBnext_combination(ar, n, k);
    }
    return result;
}

int MBnext_combination(int *ar, int n, int k) {
    int finished = 0;
    int changed = 0;
    int i, j;

    if (k > 0) {
        for (i = k - 1; !finished && !changed; i--) {
            if (ar[i] < (n - 1) - (k - 1) + i) {
                /* Increment this element */
                ar[i]++;
                if (i < k - 1) {
                    /* Turn the elements after it into a linear sequence */
                    for (j = i + 1; j < k; j++) {
                        ar[j] = ar[j - 1] + 1;
                    }
                }
                changed = 1;
            }
            finished = i == 0;
        }
        if (!changed) {
            /* Reset to first combination */
            for (i = 0; i < k; i++) {
                ar[i] = i;
            }
        }
    }
    return changed;
}

