# Challenge2 #
this is the README file for the 2nd challege of PACS corse
in the src folder there are the 5 files containing Declarations,implementation and test for the sparse Matrix. There is a file .mtx that is used for tested the class

# How run the code
In the folder src there is a MakeFile that use a the macro PACS_ROOT.
1. In case you do not have that variables before running the code on terminale you must digit export PACS_ROOT=complete_path_to_pacs-exampples_folder/Examples. for having the full path it is enought to go in our local repository pacs-examples/Examples and on terminal digits pwd
2. In case your device is a MAC book and you use docker for having the linux machine you also must set the variable LD_LIBRARY_PATH. On the terminal you have to write export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/full-path/to/the/library in our case will be the sub folder lib of pacs-examples/Examples
3. Now in Challenge2/src you can do make on terminal 

