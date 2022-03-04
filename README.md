# Portfolio
This repository compiles a few programs I've created in both C++ and Python. These projects range from physics research using C++ and the ROOT data analysis framework to simple projects done in python. The purpose of this respository is to display the skills I have developed in multiple languages. 
____________________________________________________________________________________________________________________________________________________________________
About Each Program: 

  spatial.C:
  This was one of the first programs I created as an undergraduate researcher. The program was created with C++ and is run in the Root data analysis framework on     linux. The program reads data in from files containing information from the MINOS neutrino experiment. From there, I anylize two conditions; First I calculate the 
  angular difference between neutrinos in a cosmic ray event where the number of tracks is greater than 2. Next, a similar calculation is made for the number of       tracks in an event equal to 2. The result of running this program is two histograms, one for ntracks > 2 and the other for ntracks == 2. 
  
  InterpLab:
  This project was completed as a lab for a computational physics course. The folder contains 3 files, Interpolation.py is the main program, Interpf.py is the         functions program, and Interp.pdf is the result of running the program. 
  The program is relatively simple, it finds the interpolation of a function and creates two   plots; One being a 4th degree interpolation and the other an 8th         degree interpolation.
  
  2D_Ising.py:
  This program was created for my final project in a computational physics course. The program uses the background physics of the ising model to create a simulation   using a monte carlo method supported by the metropolis algorithim. The purpose of the program is to reproduce the results of the magnetization vs. Temperature       relationship defined by the Ising model. 
  
  Bool.C:
  This is an older version of my undergraduate research project. I will eventually update this with the final version once it is completed. This program is written   in C++ and is run in the ROOT data analysis framework. The purpose of this program is to simulate incoming neutrinos on a neutrino detector. The detector is         modeled off of the MINOS far detector, it has an octagonal prism shape that is created in the prgram using plane definitions. A monte carlo method is utilized;   This  portion generates random points on a sphere each with directional vectors. The program then uses the parametric equation to determine whether the trajecotry of these randomly generated rays will hit the detector or not. The final version of this project finds the surface area of the detector as a function of angles  theta  and phi of the incoming ray.
  
  

