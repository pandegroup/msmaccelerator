Running Workers on Clusters
---------------------------

TJL 6.9.12


This is a short guide about how to set up a cluster to run workers. We've
provided some scripts to get you going immediately, though you may need
to modify these for your own personal configuration.

To get going, run

    $ autoinstall_linux.sh

which should work on almost any Linux system. This is a simple script
that downloads and installs

+ FFTW    (v 3.3.1)
+ GROMACS (v 4.5.3)
+ CCTools (v 3.4.2)

and thus gives you everything you need to get going with simple GROMACS
runs. If you want to run another MD package, you'll need to install it
yourself.

NOTE: Be a little careful using this script. Because of the vast
heterogeneity in clusters, something is likely to break, and you'll have
to make small modifications.

Once all the necessary software is installed, it's time to boot up some
workers! Submit a bunch of workers with the scripts

  sge_submit_workers
  pbs_submit_workers
  pbs_submit_workers_tunnel (see notes below)

each of which should be in your path. These scripts will each submit a
specified number of jobs, each with a number of workers running a number
of threads. E.G.

  Usage: sge_submit_workers [options] <servername> <port> <num-workers>
  $ ./pbs_submit_workers mycomputer 5521 12

submits 12 PBS jobs, each running X workers, each worker running X threads.

That's pretty much it! Boot up a master and get crackin'!

Notes on the script "pbs_submit_workers_tunnel":
This script is designed to work similarly to the others, but on clusters
where compute nodes cannot connect to the internet (but the head node can).
In this script, the worker jobs starts an ssh tunnel that routes communication
to the server/master through the cluster head node. If you're going to use
this script we strongly recommend you ensure everything is working for your
particular system.

