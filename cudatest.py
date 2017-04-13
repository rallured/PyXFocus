import traces.PyTrace as PT
import time

tstart=time.time()
PT.circularbeam(200.,512*512)
print time.time()-tstart
gpus = PT.transferToGPU()
print time.time()-tstart
PT.radgratC[512,512](12.e3,160./12.e3,1,1.,*gpus[:6])
PT.radgratC[512,512](12.e3,160./12.e3,1,1.,*gpus[:6])
PT.radgratC[512,512](12.e3,160./12.e3,1,1.,*gpus[:6])
PT.radgratC[512,512](12.e3,160./12.e3,1,1.,*gpus[:6])
PT.returnFromGPU(*gpus)
print time.time()-tstart

##PT.circularbeam(200.,512*512)
##print time.time()-tstart
##PT.radgrat(12.e3,160./12.e3,1,1.)
##PT.radgrat(12.e3,160./12.e3,1,1.)
##PT.radgrat(12.e3,160./12.e3,1,1.)
##PT.radgrat(12.e3,160./12.e3,1,1.)
##print time.time()-tstart
