/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/* 
 * @file   FunctionMap.cpp
 * @author Francis Ricci
 * 
 * Created on November 24, 2014, 4:23 PM
 */

#include "FunctionMap.h"

#include <stdexcept>

int funcs_to_am(int nFuncs){
	switch(nFuncs){
	case 1: return 0;
	case 3: return 1;
	case 5: return 2;
	case 6: return 2;
	case 7: return 3;
	case 10: return 3;
	default: throw std::invalid_argument("Invalid number of angular momentum functions used.");
	}
}

int get_am(int i, std::vector<int>& nFuncs){
	int numFuncs = 0;
	for (int i : nFuncs){
		numFuncs += nFuncs.at(i);
		if (numFuncs >= i) return funcs_to_am(i);
	}
	 throw std::invalid_argument("Invalid number of angular momentum functions used.");
}

int diag_address(int i, int j){

	if (i < j){
		int tmp = i;
		i = j;
		j = tmp;
	}

	return (i + 1) * i / 2 + j;
}

int bas_to_shell(int i, std::vector<int>& nFunc){
	int shell = 0;
	while (i > 0)
		i -= nFunc[shell++];
	if (i != 0)
		shell--;
	return shell;
}

int bas_in_shell(int i, std::vector<int>& nFunc){
	int shell = 0;
	while (i > 0){
		i -= nFunc[shell++];
	}
	if (i == 0)
		return 0;
	return i + nFunc[shell - 1];
}
