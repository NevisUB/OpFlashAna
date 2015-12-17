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
from ROOT import signalana

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify analysis output root file name
my_proc.set_ana_output_file("michel_mc.root");

# Specify data output root file name
my_proc.set_output_file("")

# muon PE thershold
muonPE = 50

# signal processing module instance
signalProcessor = signalana.SignalProcessing()
signalProcessor.setHitPEDifferentialThresh(1)#PE
signalProcessor.setDeadTime(3)#usec
signalProcessor.setMuonPEThresh(muonPE)#usec

mergewfs = fmwk.OverlayWF_Paddles()
mergewfs.setPMTProducer("saturation")
#mergewfs.setProducer("pmtreadout")
mergewfs.setTrigProducer("triggersim")
#mergewfs.setTrigProducer("daq")
mergewfs.setSignalProcessor(signalProcessor)
mergewfs.useMC(True)
mergewfs.setMuonPEThresh(muonPE)#PE
mergewfs.setNoisePE(2)#PE
mergewfs.setRequireMuonPeak(True)
mergewfs.setMaximumMuonTime(25)
mergewfs.setMaximumMuonNumber(1)
mergewfs.setLateLightAmplitude(0.25)
mergewfs.setLateLightTimeConstant(0.85)
mergewfs.setVerbose(True)
my_proc.add_process(mergewfs)

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

