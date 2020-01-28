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

ITNode* NMC(BELLApars& b_pars, ITNode *root)
{
	std::vector<Interval> intervals;
	Interval tmp;

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
	int trace = b_pars.minOverlap; 
	for(int i = 2; i < (b_pars.upperNMC+1); i++)
	{
   		if(PyCallable_Check(pFunc))
   		{
    	    double myProb = 1 - b_pars.errorRate;
   		    pValue = Py_BuildValue("(diiid)", myProb, b_pars.kmerSize, i, trace, b_pars.minpNMC);
   		    PyErr_Print();
   		    presult = PyObject_CallObject(pFunc, pValue);
   		    PyErr_Print();
			trace = PyInt_AsLong(presult);
	
   		} else PyErr_Print();

   		tmp.upper = PyInt_AsLong(presult);
		tmp.num   = i-1;

		intervals.push_back(tmp);

		tmp.lower = tmp.upper;
		Py_DECREF(pValue);
	}

	// build interval tree
    for (int i = 0; i < intervals.size(); i++)
		root = insert(root, intervals[i]);

	// inorder(root);

   	// clean up
   	Py_DECREF(pModule);
   	Py_DECREF(pName);
   	Py_Finalize();

	return root;
}

#endif //MARKOV_H