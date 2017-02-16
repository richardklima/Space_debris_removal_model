# THIS IS OLD VERSION CONSISTENT WITH GAMES PAPER FOR NEW VERSION SEE REPOSITORY space_debris_removal
Space debris removal model based on PyKEP scientific library. Developed within ESA Ariadna study project (University of Liverpool and European Space Agency cooperation).  

Authors: Richard Klima, Daan Bloembergen, Rahul Savani, Karl Tuyls, Daniel Hennes and Dario Izzo.  
Corresponding author: Richard Klima - richard.klima(at)liverpool(dot)ac(dot)uk  

This is a raw version with minimum of comments and rather unclear structure. The code is being cleaned up and commented properly and will be updated soon.  

There is a break-up model (breakup.py) and collision model (cube_2agents.py) which is based on CUBE method for evaluating probability of collision of two objects.  

You can run the model by running the main class cube_2agents.py. To run this model successfully you will need to install PyKEP scientific library available at https://esa.github.io/pykep/.  


The initial setting is:  

time horizon - 150 years  
US removes 1 risky object every two years  
EU, China and Russia do not remove any objects  

Input (publicly available):  

SATCAT - satellite catalogue  
TLE - two-line element set database  

Output:  

There are 3 output files:  
(i)   
Number of objects every year of the simulation  
Number of catastrophic and non-catastrophic collisions  
Total number of debris  

(ii)  
Collisions of important assets for each of the players  
Risks to important assets  
Number of important assets for each of the players  

(iii)  
List of removed objects  

