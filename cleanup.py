import subprocess

subprocess.call('rm run_log.txt', shell = True)
subprocess.call('rm log', shell = True)
subprocess.call('rm log2', shell = True)
subprocess.call('rm slurm*', shell = True)
subprocess.call('rm *.pyc', shell = True)
subprocess.call('rm dsmc.ini', shell = True)
subprocess.call('rm -rf state_files', shell = True)
subprocess.call('rm -rf movie_files', shell = True)
subprocess.call('rm -rf statistics', shell = True)
subprocess.call('rm -rf plots', shell = True)
subprocess.call('rm Tocontinue', shell = True)

print "Cleanup complete"