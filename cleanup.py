import subprocess

subprocess.call('rm run_log.txt', shell = True)
subprocess.call('rm log', shell = True)
subprocess.call('rm -rf state_files', shell = True)
subprocess.call('rm -rf movie_files', shell = True)
subprocess.call('rm -rf statistics', shell = True)
subprocess.call('rm Tocontinue', shell = True)
subprocess.call('rm dsmcconfig.pyc', shell = True)

print "Cleanup complete"