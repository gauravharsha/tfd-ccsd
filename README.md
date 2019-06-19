--------------------
# Background
--------------------
The `zero-limit` branch seems to be working well, and giving the appropriate Zero Temperature limit of the amplitudes (at least for some cases - need more checks). It is, therefore, time to clean up the code and prepare a public version of the same. The Fixed Reference is obsolete, except for the purpose of the upcoming paper. Besides, there are other features that are desirable in the new code. Here we will list them out.

--------------------
## Requirements
--------------------

Overall features
1. No fixed reference needed (can perhaps be added later if deemed necessary)
2. Object Oriented structure for the evolution / driver functions.
3. Except `Main` and `IO`, keep everything in Fortran.
4. *Phase I : Use `scipy.integrate.ode' to do the integration
   *Phase II: (Later) switch almost everything to Fortran.
5. Test driven development

Features for Fortran Modules
0. Re write the beta and mu evolutions for CC with the accurate analytical expressions.
1. Keep it simple and clean -- optimize moderately (without OpenMP)
2. Avoid the numerous different do loops
3. Include custom compiler flags for better / faster performance

Features of Python Modules
1. Modular and Object Oriented

--------------------
## Important
--------------------

The division by small x, y (or the Bogoliubov parameters) in the driver terms of the Beta and Mu evolution is introducint numerical noise -- better if we get cleaner analytical expressions and code them. Use drudge to obtain them.



PRB 34 5390 (1986)
PRB 62 4927(200))
