import re

def checkinst(inst, plat):
    """Checking instruments and platforms match. Python coding of
       compatibleInstrumentAndPlatform in icedrift_instruments.c"""

    instplats = {}
    instplats['ssmi'] = [r'f(\d+)']
    instplats['ssmis'] = [r'f(\d+)']
    instplats['amsr'] = [r'.*aq.*']
    instplats['amsr2'] = [r'.*gw.*']
    instplats['ascat'] = [r'.*metop.*']
    instplats['seawinds'] = [r'qscat']
    instplats['miras'] = [r'smos']
    instplats['multi'] = [r'sel', r'oi']

    errstrs = {}
    errstrs['ssmi'] = "Not a DMSP platform ({}) for an SSM/I instrument"
    errstrs['ssmis'] = "Not a DMSP platform ({}) for an SSMIS instrument"
    errstrs['amsr'] = "Not Aqua platform ({}) for AMSR instrument"
    errstrs['amsr2'] = "Not gw? platform ({}) for AMSR2 instrument"
    errstrs['ascat'] = "Not metop? platform ({}) for ASCAT instrument"
    errstrs['seawinds'] = "Not qscat platform ({}) for SEAWINDS instrument"
    errstrs['miras'] = "Not smos platform ({}) for MIRAS instrument"
    errstrs['multi'] = "Not 'oi' or 'sel' method ({}) for MULTI instrument"

    # First check that the instrument is OK
    if not inst in instplats:
        raise ValueError("Unknown instrument ({})".format(inst))

    # Check that the platform belongs with the instrument
    match = 0
    for platexpr in instplats[inst]:
        platregex = re.compile(platexpr)
        if re.fullmatch(platregex, plat):
            match = 1
            break
    if not match:
        raise ValueError(errstrs[inst].format(plat))

    # Check DMSP number for SSM/I and SSMIS
    if (inst == 'ssmi') or (inst == 'ssmis'):
        platregex = re.compile(instplats[inst][0])
        dmsp_num = re.fullmatch(platregex, plat).group(1)
        if inst == 'ssmi':
            if int(dmsp_num) >= 16:
                raise ValueError("Platform is {} (max for SSM/I is f15)"
                                 "".format(dmsp_num))
        if inst == 'ssmis':
            if int(dmsp_num) < 16:
                raise ValueError("Platform is {} (min for SSMIS is f16)"
                                 "".format(dmsp_num))


if __name__ == '__main__':

    import argparse

    p = argparse.ArgumentParser("checkinst")
    p.add_argument('inst')
    p.add_argument('plat')
    args = p.parse_args()
    inst = args.inst
    plat = args.plat

    print("Checking instrument and platform:")
    print("Instrument: ", inst)
    print("Platform:", plat)

    checkinst(inst, plat)
