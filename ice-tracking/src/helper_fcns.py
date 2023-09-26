from datetime import datetime, timedelta

def fmtdate(dt):
    """Convert a date in datetime format to YYYYmmddHH or YYYYmmddHHMMSS
       with the time of day as 12:00"""

    if datetime.strftime(dt, "%H%M%S") == '000000':
        ymdhdate = datetime.strftime(dt + timedelta(0.5), "%Y%m%d%H")
    elif datetime.strftime(dt, "%H%M%S") == '120000':
        ymdhdate = datetime.strftime(dt, "%Y%m%d%H")
    else:
        ymdhdate = datetime.strftime(dt, "%Y%m%d%H%M%S")

    return ymdhdate


def corr_time_from_unit(time_val, time_unit):
    """Convert a NetCDF time into seconds since 1970. Based on C code
       correct_time_from_unit by T. Lavergne"""

    import re

    # Check that the units are "seconds since"
    if re.match("^seconds since", time_unit.strip()) is None:
        raise ValueError("Time unit must match 'seconds since'")

    # Take all the digits to make the YYYYMMDDHHMMSS and find the
    # difference between this reference date and 1970 to get the offset
    epochdatestr = "".join(re.findall(r'\d+', time_unit))
    epochdate = datetime.strptime(epochdatestr, "%Y%m%d%H%M%S")
    date1970 = datetime.strptime("19700101000000", "%Y%m%d%H%M%S")
    offset = (epochdate - date1970).total_seconds()

    # Add the offset to the time value and round to nearest second
    return round(time_val + offset)


def round_sig(value, sig=3):
    """Round a value to a specified number of significant figures"""
    from math import log10, floor

    return round(value, sig-int(floor(log10(abs(value))))-1)
