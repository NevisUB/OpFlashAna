import sys
from ROOT import gSystem
gSystem.Load("libOpFlashAna_OpMichel")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing OpMichel..."

sys.exit(0)

