This code accompanies publication 'Assesing the impact of imperfect adherence to artemether-lumefantrine on malaria treatment outcomes using within-host modelling' (currently in press with Nature Communications).

Output from the model is saved as text files. For single-simulation runs (in both treated and untreated models), Python scripts have been provided to produce figures of the resulting parasitaemia. For Linux/Mac users, the Makefile allows the model to be compiled, run, and (if desired) output displayed on screen using a single command. In the terminal (assuming all files are saved in the same folder), typing **make para** and hitting return will run the model for the untreated infection, while typing **make para_with_PKPD** and hitting return will run the mdoel for an infection treated with six doses of artemether-lumefantrine (AL), administered at the recommended timings. For Window users familiar with Cygwin, the makefile still works (in my experience) although, in the compilation command, you may need to specify which version of gcc you are running.

To run the model that utilising the real-world adherence data you will need to download the data from [this repository](http://actc.lshtm.ac.uk). Then remove the column headings and save as a text file. We did not use data from all 659 patients: as explained in the article, we removed some individuals who had taken too many pills in one of their doses. This was done because our model does not take toxicity effects into consideration. However, we include all the data here, so that the user can make their own decision for how to deal with this issue.

<em>Update 26th Nov. 2017:</em> I have updated the code for the treated malaria infection to make the program run more quickly. As the drug artemether is absorbed and eliminated very quickly, the ODEs can be quite stiff. I have now replaced the Runge-Kutta routine for the ODEs with a function that returns the exact analytical solution of these equations at each time point. This means that the model can safely run with a larger time step.

Here is a description of each file in the project:

* WithinHostModel.cpp (C++ code for untreated malaria infection)
* plot_parasitaemia.py (Python script to graph output from above model)
* parameters.h (Header file, storing parameters common to the C++ files)
* RK4.cpp (functions for the pharmacokinetic ODE model)
* WithinHostModel_with_PKPD.cpp (C++ code for malaria infection treated with AL)
* WithinHostModel_with_PKPD_faster.cpp (Improved C++ code for malaria infection treated with AL. This is now the recommended version to use)
* plot_parasitaemia_with_PKPD.py (Python script to graph output from above model)
* WithinHostModel_with_adherence.cpp (C++ code for malaria infection treated with AL, timings of doses determined by adherence data)
* makefile (directive which, for Linux or Mac, compiles and runs C++ model, and [if desired] runs python script which produces a visualisation of the model output)