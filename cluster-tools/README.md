Running Workers on Clusters
---------------------------

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
workers! Submit a bunch of workers with the script

    pbs_submit_workers

which should be in your path. To figure out how to use this, do the obvious

	./pbs_submit_workers -h
	
yielding
	
	Use: pbs_submit_workers [options] <servername> <port> <num-nodes> <ppn> <pbs-queue> <walltime> <workers-per-node>
	where options are:
	  -a               Enable auto mode.
	  -f               Set up SSH forwarding through the head node (if work nodes cannot connect to the internet)
	  -g               Enable GPU run (sets the CUDA device ID)
	  -s               Run as a shared worker.
	  -N <name>        Preferred master name.
	  -C <catalog>     Set catalog server to <catalog>. <catalog> format: HOSTNAME:PORT.
	  -t <time>        Abort after this amount of idle time. (default=900s)
	  -p <parameters>  SGE qsub parameters.
	  -h               Show this help message.

thus, a normal command call on one of our clusters looks like

    ./pbs_submit_workers vspm42 5521 3 24 default 24:00:00 3

which starts up PBS jobs (default queue) on 3 nodes, each running 3 workers, with each worker using 8 threads (24 ppn). These workers try to communicate with the machine `vspm42` on port 5521, and die after 24 hours.

That's pretty much it! Boot up a master and get crackin'!

Notes on the `-f` flag:
On clusters where compute nodes cannot connect to the internet (but the head node can) the worker jobs starts an ssh tunnel that routes communication
to the server/master through the cluster head node. If you're going to use
this feature we strongly recommend you ensure everything is working for your
particular system.

