MSMAccelerator
==============

Overview
--------

MSMAccelerator is an adaptive sampling engine for Markov state model powered
molecular dynamics simulations.

The basic idea is to interlace simulation and analysis together --
MSMAccelerator automates the process of building models, identifying the
undersampled regions of phase space, shooting new simulations, and then
reanalyzing the data.

The MD simulations can be done remotely with heterogeneous compute hardware.
MSMAccelerator connects to the compute nodes using WorkQueue
(http://nd.edu/~ccl/software/workqueue/), a scalable Master/Worker framework.
The only requirement is that the compute nodes be able to comunicate with the
head node via an open port. MSM construction is done on the "head" node by
the MSMAccelerator executable.

MSMAccelerator currently interfaces with Gromacs and OpenMM. Support for NAMD is
in the works.

Architecture
------------

MSMAccelerator is centered around two components. First, the WorkQueue
master/worker framework is the bridge between the command and control processes
and the worker nodes running simulation. MSMAccelerator submits "jobs" to a queue,
WorkQueue handles sending them to remote workers, running the command(s) and
then returning the MD trajectories. When trajectories are returned, MSMAccelerator
is triggered and begins data analysis.

The second major component is a MySQL database that stores (paths to) trajectories,
models, and metadata associated with the running set of simulations. The various
portions of MSMAccelerator largely stay in sync by coordinating via the database.
Using MySQL is a bit of a pain compared to SQLite, but (I think) required. The
key reason is that MSMAccelerator is a threaded application in which both threads
need to make write operations, and SQLite has no support for this type of
concurrency.


There are three main modules. First is `QMaster`, which controls the WorkQueue 
instance, adding jobs to the queue and retrieving them from the queue. The second
is `Brain`, which creates the jobs, and the third is `Builder`, which builds the
Markov state models.

In addition to these components, there is a core `Project` class that loads the
configuration file (project.yaml) and initiates the database connection. The
`sampling` module contains the adaptive sampling logic.

MSMAccelerator abstracts out the difference between MD engines by using a set
of python scripts called drivers. Each MD engine is associated with a driver
that is responsible for taking a small number of command line options (PDB file,
forcefield name, etc), and then setting up and running the simulation. Details
like the temperature, type of thermostat, etc, should be set inside the driver
file and are thus invisible to the rest of MSMAccelerator.

Dependencies
------------

MSMAccelerator is written in python2.7. We recommend the Enthought python
distribution (EPD). It depends on MSMBuilder (https://simtk.org/home/msmbuilder/).

WorkQueue and its python bindings are required, and are not always trivial to
install. WorkQueue is distributed as a part of "CCTools", and we use a forked
version of the CCTools package, which is hosted at
https://github.com/rmcgibbo/cctools-3.4.2-fork. See the install instructions in
that package for details on compiling it.

For the database component, we use MySQL, mysql-python and SQLAlchemy. Either a
dedicated MySQL server or a running instance on your local machine is required.

mysql-python does not ship by default with EPD, but can be installed with

$ easy_install mysql-python

On a fresh Ubuntu 11.10 instance, I used the following commands to get mysql setup
(with python bindings)

    sudo apt-get install mysql
    sudo apt-get install libmysqlclient-dev
    easy_install mysql-python

If you're not using EPD python, the python packages that you'll need to install
include, but may not be limited to

numpy
scipy
tables
MSMBuilder
SQLAlchemy
mysql-python

