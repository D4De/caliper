/*
****************************** Copyright & License *****************************
CALIPER v1.0 is a software framework for the reliability lifetime evaluation of Multicore architectures. Copyright (C) 2014 Politecnico di Milano.

This framework is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This framework is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details (http://www.gnu.org/licenses/).

Neither the name of Politecnico di Milano nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
********************************************************************************
*/

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>

#define EMPTY_SET "#"
#define MAKE_STRING( msg )  ( ((std::ostringstream&)((std::ostringstream() << '\x0') << msg)).str().substr(1) )

void generateTree(int max, int min, std::set<std::string>& generated);

void stringTokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters);

void loadConfigurations(const std::string& fileName, std::map<std::string, std::vector<std::vector<double> > >& configurations, int& init_num_of_cores,
        int& min_num_of_cores);
