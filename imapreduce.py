import os
import logging
from IPython.parallel import Client
from ighumanizer3.extra.svm.svm_type import *

logging.basicConfig(format=u'%(filename)s[LINE:%(lineno)d]# %(levelname)-8s [%(asctime)s]  %(message)s',
                    level=logging.DEBUG, filename="test_log.log")

logging.debug("Client creation.")
rc = Client()
dview = rc[:]
logging.debug("OS import.")
dview.execute("import os").get()


cd = dview.apply(lambda: os.path.abspath(os.curdir)).get()
logging.debug("Current dir:" + str(cd))
dview.apply(lambda: os.chdir("/home/yakovlev/devel/ig-assist-tools")).get()
cd = dview.apply(lambda: os.path.abspath(os.curdir)).get()
logging.debug("Current dir (after changing):" + str(cd))
dview.execute("from ighumanizer3.extra.svm.svm_type import *").get()

dview.push({'test': test}).get()

logging.debug("Test started!")
res = dview.map(test, range(7, 11))
print(res)