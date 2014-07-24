/*
****************************** Copyright & License *****************************
CALIPER v1.0 is a software framework for the reliability lifetime evaluation of Multicore architectures. Copyright (C) 2014 Politecnico di Milano.

This framework is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This framework is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details (http://www.gnu.org/licenses/).

Neither the name of Politecnico di Milano nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
********************************************************************************

 This code generates the aging rate configuration file given the topology of a multicore system, based on a simple thermal model
 */

#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <set>
#include "utils.h"

#define ROWS 5
#define COLS 5

// Electro-Migration related parameters
#define BETA 2
#define ACTIVATIONENERGY 0.48
#define BOLTZMANCONSTANT 8.6173324*0.00001
#define CONST_JMJCRIT 1500000
#define CONST_N 1.1
#define CONST_ERRF 0.88623
#define CONST_A0 30000 //cross section = 1um^2  material constant = 3*10^13

// Thermal model parameters
#define ENV_TEMP 295 //room temperature
#define SELF_TEMP 40 //self contribution
#define NEIGH_TEMP 5 //neighbor contribution

#define getAlpha(temp) ((CONST_A0 * (pow(CONST_JMJCRIT,(-CONST_N))) * exp(ACTIVATIONENERGY / (BOLTZMANCONSTANT * temp))) / CONST_ERRF)

// Temperature model computation
void tempModel(double loads[][COLS], double temps[][COLS]);

void tempModel(double loads[][COLS], double temps[][COLS], int rows, int cols) {
    double temp;
    int i, j, k, h;
    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++) {
            for (k = -1, temp = 0; k < 2; k++)
                for (h = -1; h < 2; h++)
                    if ((k != 0 || h != 0) && k != h && k != -h && i + k >= 0 && i + k < rows && j + h >= 0 && j + h < cols){
                        temp += loads[i + k][j + h] * NEIGH_TEMP;
                    }
            temps[i][j] = ENV_TEMP + loads[i][j] * SELF_TEMP + temp;
        }
}

int main(int argc, char* argv[]) {
    int units, min_units;
    std::set<std::string> tree;
    int rows, cols;
    double loads[ROWS][COLS];
    double temps[ROWS][COLS];
    double wl;
    int QoSnotOkArch = 0;


    ////////////////////////////////////////////////////////////////////////////////
    //parsing input arguments
    ////////////////////////////////////////////////////////////////////////////////

    if (argc != 6) {
        std::cerr << "USAGE: " << argv[0] << " <rows> <cols> <min_unit> <initialWorkloadXunit> <output_file>" << std::endl;
        return 0;
    }
    if (atoi(argv[1]) < 1 || atoi(argv[1]) > ROWS) {
        std::cerr << "Not valid number of rows: " << argv[1] << std::endl;
        return 0;
    }
    if (atoi(argv[2]) < 1 || atoi(argv[2]) > COLS) {
        std::cerr << "Not valid number of cols: " << argv[2] << std::endl;
        return 0;
    }
    if (atoi(argv[3]) <= 0 || atoi(argv[3]) > atoi(argv[1]) * atoi(argv[2])) {
        std::cerr << "Not valid minimum number of units: " << argv[3] << std::endl;
        return 0;
    }
    if (atof(argv[4]) < 0) {
        std::cerr << "Not valid initial workload per unit: " << argv[4] << std::endl;
        return 0;
    } else if (atof(argv[4]) > 1)
        std::cerr << "Excessive initial workload per unit: " << argv[4] << std::endl;

    char* filename = argv[5];


    ////////////////////////////////////////////////////////////////////////////////
    //set up environment
    ////////////////////////////////////////////////////////////////////////////////

    rows = atoi(argv[1]);
    cols = atoi(argv[2]);
    wl = atof(argv[4]);
    units = rows * cols;
    min_units = atoi(argv[3]);

    generateTree(units, min_units, tree);

    std::cout << std::setprecision(10);
    std::cerr << std::setprecision(10);

    ////////////////////////////////////////////////////////////////////////////////
    //compute aging rates
    ////////////////////////////////////////////////////////////////////////////////
    std::ofstream outputFile(filename);
    if (outputFile.is_open()){
        for (std::set<std::string>::iterator setIt = tree.begin(); setIt != tree.end(); setIt++) {
            std::vector<std::string> tokens;
            std::set<int> failed;
            stringTokenize(*setIt, tokens, ",");
            if (tokens.at(0) != EMPTY_SET) {
                for (int i = 0; i < tokens.size(); i++) {
                    failed.insert(atoi(tokens[i].c_str()));
                }
            }
            double distributedLoad = wl * units / (units - failed.size());
            if (distributedLoad > 1 && QoSnotOkArch < units - failed.size()) {
                QoSnotOkArch = units - failed.size();
            }
            //std::cerr << *setIt << "-" << wl << "-" << distributedLoad;
            for (int i = 0, k = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    if (failed.count(k) == 0) {
                        loads[i][j] = distributedLoad;
                    } else {
                        loads[i][j] = 0;
                    }
                    //std::cerr << " " << loads[i][j];
                    k++;
                }
            }
            //std::cerr << std::endl;
            tempModel(loads, temps, rows, cols);
            outputFile << *setIt;
            double max = 0;
            for (int i = 0, k = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    if (failed.count(k) == 0)
                        outputFile << " " << getAlpha(temps[i][j]);
                    else
                        outputFile << " 0";
                    k++;
                    //if(failed.count(k)==0 && temps[i][j] > max)
                    //  max = temps[i][j];
                    //k++;
                }
            }

            //for(int i=0, k=0; i<rows; i++){
            //  for(int j=0; j<cols; j++){
            //    std::cout << " " << getAlpha(max);
            //  }
            //}
            outputFile << std::endl;        
        }
        outputFile.close();
    } else {
        std::cerr << "Error while opening the output file: " << filename << std::endl;
        exit(1);
    }

    if (QoSnotOkArch > 0)
        std::cerr << "QoS not satisfied with less than " << (QoSnotOkArch + 1) << " cores" << std::endl;

    return 0;
}

