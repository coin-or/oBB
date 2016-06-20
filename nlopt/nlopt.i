// -*- C++ -*-

%define DOCSTRING
"NLopt is a multi-language library for nonlinear optimization (local or
global, with or without derivatives, and supporting nonlinear
constraints).  Complete documentation, including a Python tutorial,
can be found at the NLopt web page: http://ab-initio.mit.edu/nlopt"
%enddef

%module(docstring=DOCSTRING) nlopt
%{
#include "nlopt.hpp"
%}

%include "std_vector.i"
namespace std {
  %template(nlopt_doublevector) vector<double>;
};

%include "nlopt-exceptions.i"

%include "nlopt-python.i"

%include "nlopt.hpp"
