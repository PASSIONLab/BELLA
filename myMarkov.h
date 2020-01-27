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
    int lower, upper;
}; 
  
// node struct in interval search tree 
struct ITNode 
{
    Interval *i;
    int max, num; 
    ITNode *left, *right; 
}; 
  
// function to create a new interval search tree node 
ITNode * create(Interval i, int num) 
{ 
    ITNode *temp = new ITNode; 
    temp->i = new Interval(i); 
    temp->max = i.upper;
	temp->num = num;
    temp->left = temp->right = NULL; 
}; 

// function to insert a new interval search tree node 
ITNode *insert(ITNode *root, Interval i, int num) 
{ 
    // tree is empty, new node becomes root 
    if (root == NULL) 
        return create(i, num); 
  
    // low value of interval at root 
    int l = root->i->lower; 
  
    if (i.lower < l) root->left = insert(root->left, i, num);
    else root->right = insert(root->right, i, num);
  
    return root; 
}

// function to check if L is in overlap
bool within(Interval i, int L) 
{ 
    if (i.lower <= L && L <= i.upper) return true; 
    return false; 
} 
  
// main function that searches a given interval i in a given interval tree 
Interval *search(ITNode *root, int L) 
{ 
    // tree is empty 
    if (root == NULL) return NULL; 
  
    if (within(*(root->i), L)) 
        return root->i->num; 
      
    if (root->left != NULL && L < root->left) 
        return search(root->left, L); 
  
    // else interval can only overlap with right subtree 
    return search(root->right, L); 
} 

int myMarkovFunc(BELLApars& b_pars)
{
	// https://stackoverflower.com/questions/3286448/calling-a-python-method-from-c-c-and-extracting-its-return-value/24687260	
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
   	pFunc = PyDict_GetItemString(pDict, (char*)"myOverlap");

   	if(PyCallable_Check(pFunc))
   	{
        double myProb = 1-b_pars.errorRate;
   	    pValue = Py_BuildValue("(di)", myProb, b_pars.kmerSize);
   	    PyErr_Print();
   	    presult = PyObject_CallObject(pFunc, pValue);
   	    PyErr_Print();
    
   	} else 
   	{
   	    PyErr_Print();
   	}

   	int result = PyInt_AsLong(presult);
   	Py_DECREF(pValue);

   	// clean up
   	Py_DECREF(pModule);
   	Py_DECREF(pName);
   	Py_Finalize();
    
    return result;
}

#endif //MARKOV_H