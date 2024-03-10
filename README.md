# Install HOOMD
The HOOMD version used here is [**3.5.0**](https://hoomd-blue.readthedocs.io/en/v3.5.0/). 

# Get parameter file
If parameter file is not avalable get it by going into param folder and then
# Set up parameter
The parameter space is included in "parentparam.yml", the parameter could be either ON/OFF, or string or a list of values.

# Generate folder with given parameters
Run 
$ py param.py parentparam.yml.

It will iterate all parameters and generate folders with json file included.

# Run simulation
$./runscript.sh all/allnew/"folder name"

Options: 
- all: run all folders.
- allnew: run all new generated folders which contains only json file.
- "folername": run given folder.

If permission is required, run

$ chmod +x runscript.sh

before run the .sh file.
