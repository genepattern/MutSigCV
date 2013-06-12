# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#

import sys
import subprocess
import os
import shutil
import glob

if len(sys.argv) < 2:
    print "Usage:"
    print "   run_matlab.py [--with-display] <libdir> <compiledMFile> args"
    sys.exit(1)

argIndex = 1
xvncProcess = None
try:
    if sys.argv[argIndex] == "--with-display":
        argIndex = argIndex + 1
        jobId = os.environ.get("LSB_JOBID")
        display = "1"
        if jobId is not None:
            display = repr((int(jobId) % 10000) + 10000)
        xvncProcess = subprocess.Popen(["Xvnc", ":" + display, "-depth", "16"])
        os.environ["DISPLAY"] = "localhost:" + display
    
    exe_dir = sys.argv[argIndex] 
    print "libdir:", exe_dir
    argIndex = argIndex + 1
    
    
    matlab=subprocess.Popen(["which", "matlab"], stdout=subprocess.PIPE).communicate()[0]
    
    matlabHome=os.path.dirname(os.path.dirname(matlab))
    LD_LIBRARY_PATH=".:%s/runtime/glnxa64" % (matlabHome)
    LD_LIBRARY_PATH="%s:%s/bin/glnxa64" % (LD_LIBRARY_PATH, matlabHome)
    LD_LIBRARY_PATH="%s:%s/sys/os/glnxa64" % (LD_LIBRARY_PATH, matlabHome)
    MCRJRE=matlabHome + "/sys/java/jre/glnxa64/jre/lib/amd64"
    LD_LIBRARY_PATH="%s:%s/native_threads" % (LD_LIBRARY_PATH, MCRJRE) 
    LD_LIBRARY_PATH="%s:%s/server" % (LD_LIBRARY_PATH, MCRJRE)
    LD_LIBRARY_PATH="%s:%s/client" % (LD_LIBRARY_PATH, MCRJRE)
    LD_LIBRARY_PATH="%s:%s" % (LD_LIBRARY_PATH, MCRJRE)
    
    XAPPLRESDIR=matlabHome + "/X11/app-defaults"
    os.environ["LD_LIBRARY_PATH"] = LD_LIBRARY_PATH
    os.environ["XAPPLRESDIR"] = XAPPLRESDIR
    print "LD_LIBRARY_PATH is", LD_LIBRARY_PATH
    
    MCR_CACHE_ROOT = os.getcwd()
    os.environ["MCR_CACHE_ROOT"] = MCR_CACHE_ROOT
    print "MCR_CACHE_ROOT is ", MCR_CACHE_ROOT
    
    exe_name = sys.argv[argIndex]
    print "mFile:", exe_name
    argIndex = argIndex + 1
    
    exe_path = os.path.join(exe_dir, exe_name)
    subprocess.check_call(["chmod", "+x", exe_path])
    os.environ["PATH"] = ".:" + os.environ["PATH"]
    matlabCommand = [exe_path] + sys.argv[argIndex:]
    print "executing:", " ".join(matlabCommand)
    subprocess.check_call(matlabCommand)
finally:
    if xvncProcess is not None:
        xvncProcess.terminate()
    mcrCacheName = ".mcrCache*"
    mcrCacheFiles = glob.glob(mcrCacheName)
    for mcrCache in mcrCacheFiles:
        if os.path.isdir(mcrCache):
            shutil.rmtree(mcrCache)
        else:
            os.remove(mcrCache)