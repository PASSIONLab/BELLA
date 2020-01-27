#ifndef MARKOV_H
#define MARKOV_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <string>
#include <stdlib.h>
#include <utility>
#include <Python.h>
#include <cstring>
#include <string.h>
#include <cassert>
#include <ios>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/sysctl.h>

#include "mtspgemm2017/common.h"

// modified from https://www.geeksforgeeks.org/interval-tree/
struct Interval 
{ 
    int lower, upper, num;
}; 
  
// node struct in interval search tree 
struct ITNode 
{
    Interval *i;
    ITNode *le, *ri; 
}; 
  
// function to create a new interval search tree node 
ITNode * create(Interval i) 
{ 
	ITNode *temp = new ITNode; 
	temp->i = new Interval(i);
	temp->le = temp->ri = NULL; 
}; 

// function to insert a new interval search tree node 
ITNode *insert(ITNode *root, Interval i) 
{ 
    // tree is empty, new node becomes root 
    if (root == NULL) 
        return create(i); 
  
    // low value of interval at root 
    int l = root->i->lower; 
  
    if (i.lower < l) root->le = insert(root->le, i);
    else root->ri = insert(root->ri, i);
  
    return root; 
}

// function to check if L is in overlap
bool within(Interval i, int L) 
{ 
    if (i.lower <= L && L < i.upper) return true; 
    return false; 
}

// main function that searches a given interval i in a given interval tree 
Interval *search(ITNode *root, int L) 
{ 
    // Base case the tree is empty
    if (root == NULL) return NULL; 
  
	// L is in this interval
    if (within(*(root->i), L)) 
		return root->i; 

	// L is greater and we are in the ri-most node
	if (root->ri == NULL && L >= root->i.upper) 
	{
		root->i.num++;
		return root->i;	
	}

	// L is smaller and we are can go left
    if (root->le != NULL && L < root->i.lower) 
		return search(root->le, L);
	
	// L is greater and we are can go left
	if (root->ri != NULL && L >= root->i.upper) 
		return search(root->ri, L); 
}

void NMC(BELLApars& b_pars, ITNode* root, int n)
{
	std::vector<Interval> intervals;
	Interval tmp;
    root = NULL; 

	// https://stackoverflow.com/questions/3286448/calling-a-python-method-from-c-c-and-extracting-its-return-value/24687260	
  	// set PYTHONPATH to working directory
   	setenv("PYTHONPATH",".",1);
   	PyObject *pName, *pModule, *pDict, *pFunc, *pValue, *presult;
   	// init
   	Py_Initialize();
   	// build the name object
   	pName = PyString_FromString((char*)"myMarkov");
   	// load the module object
   	pModule = PyImport_Import(pName);
   	// pDict is a borrowed reference 
   	pDict = PyModule_GetDict(pModule);
   	// pFunc is also a borrowed reference 
   	pFunc = PyDict_GetItemString(pDict, (char*)"boundary");

	tmp.lower = 0;
	for(int i = 2; i < n+1; i++)
	{
   		if(PyCallable_Check(pFunc))
   		{
    	    double myProb = 1 - b_pars.errorRate;
   		    pValue = Py_BuildValue("(di)", myProb, b_pars.kmerSize, i, b_pars.minOverlap, b_pars.minpNMC);
   		    PyErr_Print();
   		    presult = PyObject_CallObject(pFunc, pValue);
   		    PyErr_Print();
	
   		} else PyErr_Print();

   		tmp.upper = PyInt_AsLong(presult);
		tmp.num   = i-1;
		intervals.push_back(tmp);

		tmp.lower = tmp.upper;
		Py_DECREF(pValue);
	}

   	// clean up
   	Py_DECREF(pModule);
   	Py_DECREF(pName);
   	Py_Finalize();

	// build interval tree
    for (int i = 0; i < n-1; i++)
        root = insert(root, intervals[i]);
}

#endif //MARKOV_H