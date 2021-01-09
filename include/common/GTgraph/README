/* 
 * GTgraph: A suite of synthetic graph generators
 * Copyright (C) 2006  Kamesh Madduri, David A. Bader 
 * 
 * Last Modified: Feb 19, 2006
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
 * USA.
 *
 * Authors:         
 *			David A. Bader   <http://www.cc.gatech.edu/~bader>
 *                      College of Computing
 *                      Georgia Institute of Technology
 *                      Atlanta, GA 30332
 *
 *                      Kamesh Madduri   <http://www.cc.gatech.edu/~kamesh>
 *			College of Computing
 *                      Georgia Institute of Technology
 *                      Atlanta, GA 30332
 *
 */


 *******
  About 
 *******
  
  This package contains three synthetic graph generators:
  
  a) SSCA2: This generator produces graphs used in the DARPA HPCS 
     SSCA2 benchmark. The graphs are directed with integer weight edges.
     For a detailed description of the graph parameters, refer
     gen.pdf in the main directory
 
  b) random: We include two random graph generators --
     i) an Erdos Renyi graph generator which takes as input arguments 
     the number of vertices (n), and the constant probability (p) of an 
     edge between any pair of vertices in the graph.  
     ii) a random graph generator which takes as input the number of 
     vertices (n) and edges (m), and adds m edges randomly choosing a 
     pair of vertices each time.

  c) R-MAT: To generate graphs with power-law degree distributions
     and small-world characteristics, we apply the Recursive Matrix
     (R-MAT) graph model discussed in [1].  

 *****************
  Getting Started		 
 *****************

 a) Set the following variables in Makefile.var
    i.  CC      - the C compiler
    ii. Verify if MAKE, CFLAGS, LDFLAGS are specified correctly
 
 b) This package includes a stripped-down version of the SPRNG2.0a
    library from Mike Mascagni and Ashok Srinivasan for portable 
    random number generation. We use the combined multiple recursive 
    generator (cmrg) which is a combination of a 64-bit LCG (lcg64) and
    the 32-bit multiple recursive generator
    Edit the SPRNG-specific variables in Makefile.var
    i.   Set the path to the Fortran compiler
    ii.  Verify other variables
    
 c) Build SPRNG and all the graph generators
    % make

    To build each generator (and SPNRG) separately, do  
    % make sprng
    % make ssca
    % make rand
    % make rmat

    For cleaning up object files before new builds
    % make clean  -- removes all object files
    % make clean-sprng      
    % make clean-ssca
    % make clean-rand
    % make clean-rmat

 d) To build and run all graph generators with the default options, do
    % make test-run
    
    Note: The default parameters are set to generate "large graphs" -- 
	  graphs with vertices and edges in the order of hundreds of 
	  millions. The output edge lists are written to file in the 
	  plain text DIMACS graph format. So please make sure that you
	  have sufficient system resources before this run.

 *******
  Usage
 *******
 
  SSCA2
 *******
 
 % GTgraph-ssca2 [-options]
	-s ###  SCALE value (integer) to use (default -- 20)
        -c ###  config file to use
        -o ###  output file to write the graph to (default -- sample.gr)
	-h      display this message

  random
 ********

 % GTgraph-random [-options]
        -c ###  config file to use
                (default: read parameters from init.c)
        -t ###  random graph model
                (0 for G(n, p), 1 for G(n, m))
        -n ###  no. of vertices
        -p ###  probability of an edge (graph model 0)
        -m ###  no. of edges (graph model 1)
        -o ###  output file to write the graph to
        -h      display this message

  R-MAT
 *******

 % GTgraph-rmat [-options]
	-c  ###	 config file to use
        -n  ###  no. of vertices (default -- 10^7)
	-m  ###	 no. of edges (default -- 10^8)
	-o  ###  output file to write the graph to
	-h       display this message

 Use the config file option to specify custom values of a, b, c and d, the 
 R-MAT generator parameters [1]. 

 For all the generators, if the executable is called 
 without any arguments, the default input parameters 
 specified in init.c are assumed.

   References
  ************
  
  [1] D. Chakrabarti, Y. Zhan and C. Faloutsos, R-MAT: A Recursive
      Model for Graph Mining, Proc. SIAM Intl. Conf. on Data Mining, 2004

