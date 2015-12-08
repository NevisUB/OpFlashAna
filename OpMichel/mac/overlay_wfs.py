import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


from ROOT import gSystem,TMath
from ROOT import larlite as fmwk
from ROOT import larutil

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify analysis output root file name
my_proc.set_ana_output_file("michels.root");

# Specify data output root file name
my_proc.set_output_file("")

mergewfs = fmwk.OverlayWF_Paddles()
mergewfs.setPMTProducer("saturation")
#mergewfs.setProducer("pmtreadout")
#mergewfs.setTrigProducer("triggersim")
mergewfs.setTrigProducer("daq")
mergewfs.useMC(False)
mergewfs.setMuonPEThresh(12)#PE
mergewfs.setHitPEDifferentialThresh(1)#PE
mergewfs.setBaselinePE(2)#PE
mergewfs.setDeadTime(3)#usec
mergewfs.setRequireMuonPeak(True)
mergewfs.setMaximumMuonTime(0.2)
mergewfs.setMaximumMuonNumber(100)
mergewfs.setVerbose(False)
my_proc.add_process(mergewfs)

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

