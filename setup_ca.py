import matplotlib
from distutils.core import setup
import py2exe

setup(windows=[{"script":"cluster_analysis.py"}],
      options={"py2exe": {"includes": ["sip", "matplotlib.backends",
                            "matplotlib.backends.backend_qt4agg",
                            "pylab", "numpy",
                            "matplotlib.backends.backend_tkagg",
                            "scipy.sparse.csgraph._validation"],
                          'bundle_files': 1,
                          'compressed': True}},
      data_files=matplotlib.get_py2exe_datafiles())