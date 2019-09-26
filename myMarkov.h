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

int myMarkovFunc(BELLApars& b_pars)
{
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