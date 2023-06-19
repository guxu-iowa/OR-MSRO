# OR-MSRO

## Background
The repository hosts the codes developed to implement decision rules that are used to solve multi-stage robust optimization problems. The decision rules are developed in the paper "Improved Decision Rule Approximations for Multi-Stage Robust Optimization via Copositive Programming". 

All optimization problems are solved using MOSEK 8.1.0.56 via the YALMIP interface on a 16-core 3.4 GHz Linux PC with 32 GB RAM. Computation times are reported in the paper. 

The repository has three folders with one for each of the three examples in the paper.

## Newsvendor 
File test.m is the entrance of the program. It generates test instances and run different decision rules over those instances. Decision rules include quadratic decision rule (newsvendor_qdr.m), copositive-based decision rule (newsvendor_cop.m), approximate s_lemma (newsvendor_as.m), and polynomial decision rule (newsvendor_poly.m). Auxiliary files including genBC.m, genind.m,and genpow.m are used to implement polynomial decision rule. 

## Inventory_control
File test.m is the entrance of the program. The file calls the file "run_decision_rules.m", which actually generates inventory control instances and compare different decision rules over these instances. The compared decision rules are implemented in approxSL.m (aproximate S-Lemma decision rule), ldr.m (linear decision rule), piecewiseLDR.m (piecewise linear decision rule) and pdr.m (polynomial decision rule) with degrees 2 and 3. 

## Index_tracking
File test.m is the entrance of the program. The file calls the file "run_decision_rules.m", which actually generates index tracking  instances and compare different decision rules over these instances. The compared decision rules are implemented in approxSL.m (aproximate S-Lemma decision rule), lqdr.m (linear and quadratic decision rule), and pdr.m (polynomial decision rule) with degrees 2 and 3. 


## References
ApS M (2016) The MOSEK optimization toolbox for MATLAB manual. Version 8.0. URL http://docs.
mosek.com/8.0/toolbox/index.html

Lofberg J (2004) Yalmip: A toolbox for modeling and optimization in matlab. Computer Aided Control
Systems Design, 2004 IEEE International Symposium on, 284–289 (IEEE)

