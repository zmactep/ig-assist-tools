import os
import logging
from IPython.parallel import Client
from ighumanizer3.extra.svm.svm_type import *

logging.basicConfig(format=u'%(filename)s[LINE:%(lineno)d]# %(levelname)-8s [%(asctime)s]  %(message)s',
                    level=logging.DEBUG, filename="test_log.log")

rc = Client()
dview = rc[:]
dview.execute("import os").get()

dview.apply(lambda: os.path.abspath(os.curdir)).get()
dview.apply(lambda: os.chdir("/home/yakovlev/devel/ig-assist-tools")).get()
dview.apply(lambda: os.path.abspath(os.curdir)).get()

dview.execute("from ighumanizer3.extra.svm.svm_type import *").get()

dview.push({'test': test}).get()

res = dview.map(test, range(7, 11))
