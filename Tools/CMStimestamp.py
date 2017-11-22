###
# Reading CMS TimeStamp
#######

fname = "/eos/cms/store/group/dpg_dt/comm_dt/dtRootple2016/Run2016HZMuPromptReco-v2.root"

file = TFile(fname,"o")

tree = file.Get("DTTree")
#tree.Print()

#  |             32 Bit                |       32 bit      |
#  |--- Time in seconds UNIX format ---|--- millisecond ---|
#

# Use the mask to obtain the millisecond part
mask = 0xFFFFFFFF

tree.GetEntry(1)
time = int(tree.timestamp)
# millisecond part
time_micro  = time & mask
# Shift by 32bit to obtain the time in seconds
time_second = time >> 32

# printing
print 'Binary: {:019b} - Decimal: {:19d} -- Decimal after 32bit shift:  {:19d}'.format(time,time,time_second)

# Conversion in date
print(
    datetime.datetime.fromtimestamp(
        int(time_second)
    ).strftime('%Y-%m-%d %H:%M:%S')
)
