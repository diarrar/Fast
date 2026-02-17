#=============================================================================
# FastS requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore
#=============================================================================
import os
from setuptools import setup, Extension
from importlib.util import spec_from_file_location, module_from_spec

def loadModuleFromPath(modname):
    # Load a Python file by filesystem path (PEP-517 isolated build requirement)
    helper = os.path.join(os.path.dirname(__file__), modname + ".py")
    spec = spec_from_file_location(modname, helper)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

# Compiler settings must be set in installBase.py / installBaseUser.py
Dist = loadModuleFromPath('../../../Cassiopee/Cassiopee/KCore/Dist')
installBase = loadModuleFromPath('../../../Cassiopee/Cassiopee/KCore/installBase')
Dist.setConfigDict(installBase.installDict)
additionalLibPaths = Dist.getAdditionalLibPaths()
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Test if xcore exists =======================================================
(xcoreVersion, xcoreIncDir, xcoreLibDir) = Dist.checkModuleCassiopee("XCore")

# Test if connector exists =====================================================
(connectorVersion, connectorIncDir, connectorLibDir) = Dist.checkModuleCassiopee("Connector")

# Test if fast exists =======================================================
#(fastcVersion, fastcIncDir, fastcLibDir) = Dist.checkFastC()
(fastcVersion, fastcIncDir, fastcLibDir) = Dist.checkModuleFast("FastC")

# Test if fasts exists =====================================================
#(fastsVersion, fastsIncDir, fastsLibDir) = Dist.checkFastS()
(fastsVersion, fastsIncDir, fastsLibDir) = Dist.checkModuleFast("FastS")

# Test if fastp exists =====================================================
#(fastpVersion, fastpIncDir, fastpLibDir) = Dist.checkFastP()

# Test if fastlbm exists =====================================================
#(fastlbmVersion, fastlbmIncDir, fastlbmLibDir) = Dist.checkFastLBM()
#(fastaslbmVersion, fastaslbmIncDir, fastaslbmLibDir) = Dist.checkFastASLBM()


# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi(additionalLibPaths, additionalIncludePaths)

# Compilation des fortrans ====================================================
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir, xcoreLibDir, connectorLibDir, fastcLibDir, fastsLibDir]
includeDirs = [numpyIncDir, kcoreIncDir, xcoreIncDir, connectorIncDir, fastcIncDir, fastsIncDir]
libraries = ["fast", "fasts", "fastc", "connector", "xcore", "kcore"]
#libraryDirs += [fastaslbmLibDir]
#includeDirs += [fastaslbmIncDir]
#libraries += ["fastaslbm"]

#libraryDirs = ["build/"+prod, kcoreLibDir, xcoreLibDir, connectorLibDir, fastcLibDir,fastsLibDir, fastlbmLibDir]
#includeDirs = [numpyIncDir, kcoreIncDir, xcoreIncDir, connectorIncDir, fastcIncDir, fastsIncDir, fastlbmIncDir]
#libraries = ["fast", "fasts", "fastlbm", "fastc", "connector", "xcore", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs

ADDITIONALCPPFLAGS=[]
if mpi:
    libraryDirs.append(mpiLibDir)
    includeDirs.append(mpiIncDir)
    ADDITIONALCPPFLAGS = ['-D_MPI']
    libraries += mpiLibs

# Extensions ==================================================================
listExtensions = []
listExtensions.append(
    Extension('Fast.fast',
              sources=['Fast/fast.cpp'],
              include_dirs=[".","Fast"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Fast",
    version="4.0",
    description="Fast NS and LBM solvers.",
    author="ONERA",
    url="https://fast.onera.fr",
    packages=['Fast'],
    package_dir={"":"."},
    ext_modules=listExtensions
)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
