# EU Surface Temperature for All Corners of Earth

**Horizon 2020**
**H2020-EO-2014 New ideas for Earth-relevant Space Applications**

EUSTACE: EU Surface Temperature for All Corners of Earth https://www.eustaceproject.org/

![EUSTACE Logo](/images/eustace.jpg)

The EUSTACE project produced:

* Surface air temperature estimates (with estimates of uncertainty) for all surfaces of Earth, derived from satellite surface skin temperature retrievals
* Global daily analyses of surface air temperature (with estimates of uncertainty) since 1850, based on combined information from satellite and in situ data sources

The data products are available via the [CEDA archive](https://catalogue.ceda.ac.uk/uuid/a52b2cc065a847b8a77a93896880349f).

The code used to produce these data products is located in the [eustace-data/eustace-system](https://github.com/eustace-data/eustace-system) repository.

# eustace-ambitious

This repository contains code modules for a partial implementation of an extended method for global daily analyses. It its present form it does not produce output.

This repo is archived and static - we are not able to reply to comments, suggestions for improvements, bug fixes or pull requests.

Further development of this code will take place under https://github.com/finnlindgren/

# Build instructions

The code uses the cmake build system:

* mkdir build
* cd build
* cmake ../src
* make

Partial system requirements under Ubuntu:
cmake,
libcppunit,
libopenmpi-dev,
libboost-all-dev,
libpnetcdf-dev,
doxygen

# Code structure

The are three main modules:

* eustace: reading the intermediate raw binary data files produced by eustace-system
* fmesher: handling unstructured triangulation mesh data structures; Originated by and
           adapted by [Finn Lindgren](https://github.com/finnlindgren/), originally part
	   of [R-INLA](https://bitbucket.org/hrue/r-inla/src/default/fmesher/)
* ambitious: A collection of submodules for an extended method, handling join estimation
     of daily mean temperatures and daily temperature ranges with non-Gaussian marginal distributions.

The ambitious C++ modules are:

* _qtool_ Matrix methods for block calculations; solving linear systems and computing conditional covariances.
* _poq_ Functions for the POwer Quantile model used for the diurnal temperature range distributions.
* _po_tool_ Program Options methods to allow algorithm and data locations options to be set in text configuration files, based on _Boost::program_options_.
* _bidirmap_ Bidirectional maps; low level methods needed e.g. to navigate the hierarchically structured global triangulation at different subdivision levels.
* _timer_ Timing tool used for algorithm testing.
* _versions_ Method for keeping track of the source file versions.
* _test_ A collection of unit tests.
* _rarandom_ Repeatable random access random number sequences: allows rerunning the exact same simulation generation, including the problem setup code, for the entire ensemble.
* _po_ The Program Options methods that are specific to the Ambitious system.
* _observe_ Methods for interfacing with the EUSTACE raw binary data reader, restructuring the data into the required internal formats.
* _debuglog_ Methods for logging program status and actions.
* _blocks_ The core system for scanning through data and structuring the model into its constituent components.
* _to_string_ Helper methods for I/O.
* _ostreamer_ Helper methods for I/O.
* _ncdf_ Wrapper module for calling plain netcdf reader functions, or handling parallel netcdf actions.
* _meshes_ Methods for interfacing with the point locator methods in _fmesher_, and wrapping the logic for space-time basis functions of various types.
* _informal_ Standalone programs that make use of the other modules, used for high-level testing, and partial implementation of the block-scale-solver needed by the PCG iterations.


EUSTACE has received funding from the European Union's Horizon 2020 Programme for Research and Innovation, under Grant Agreement no. 640171

![EC Logo](/images/logo-ce-horizontal-en-quadri-hr.jpg)
