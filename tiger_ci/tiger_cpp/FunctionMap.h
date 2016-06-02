/* 
 * @file   FunctionMap.h
 * @author Francis Ricci
 *
 * Created on November 24, 2014, 4:23 PM
 */

#ifndef FUNCTIONMAP_H
#define	FUNCTIONMAP_H

#include <vector>

/* @brief convert the number of basis functions to angular momentum
 * @param nFuncs   number of functions in this angular momentum component
 * @return size_t  angular momentum component
 */
int funcs_to_am(int nFuncs);

/* @brief find the angular momentum of a given basis function
 * @param i       the basis function index
 * @param nFuncs  basis functions per shell
 * @return size_t angular momentum component
 */
int get_am(int i, std::vector<int>& nFuncs);

/* @brief convert a set of two indices to a diagonal address
 * @param i     index 1
 * @param j     index 2
 * @return int  diagonal index
 */
int diag_address(int i, int j);

/* @brief convert a basis function index to shell index
 * @param i     basis function index
 * @param nFunc number of basis functions per shell
 * @return int  shell index
 */
int bas_to_shell(int i, std::vector<int>& nFunc);

/* @brief convert a basis function index to its index in shell
 * @param i     basis function index
 * @param nFunc number of basis functions per shell
 * @return int  index in shell
 */
int bas_in_shell(int i, std::vector<int>& nFunc);

#endif	/* FUNCTIONMAP_H */

