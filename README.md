# LDR Solver

This repository contains the active set method used in the numerical experiments in the paper:

	H. Lu and B. Sturt (2023). On the Sparsity of Optimal Linear Decision Rules in Robust Optimization


## Citation

If you use the code in this repository in your own research, please cite the above paper as follows:

```
@article{lusturt2023sparsity,
  title={On the Sparsity of Optimal Linear Decision Rules in Robust Inventory Management},
  author={Lu, Haihao and Sturt, Bradley},
  year={2023}
}
```


## Organization

This repository is organized into the following directories: 

* `examples/` Contains the scripts for performing the numerical experiments with active set method. The output of executing the scripts is saved in the `output/` directory. 
	*  `bental.jl`  After loading the script into the Julia terminal, the code is executed using the main function. 
*  `src/` Contains the code for the acive set method. 
	*  `algorithm_full.jl` Contains the code for the active set method.
  	*  `lp_solvers.jl` Contains the code for solving the robust counterpart
  	*  `utils.jl` Contains the code for detecting nonzero coefficients in the robust optimization problem. 

## Requirements

The scripts from the directory `julia_scripts` were tested using the Julia programming language (version 1.5.2) with the following packages and versions: 
* JuMP v1.3.1
* Gurobi v0.11.3. 

The author assumes responsibility for any errors in the code and paper.  

## License

MIT License

Copyright (c) 2023 Bradley Sturt

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


